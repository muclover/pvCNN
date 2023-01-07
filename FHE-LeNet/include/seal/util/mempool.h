#pragma once

#include <cstdint>
#include <cstring>
#include <vector>
#include <stdexcept>
#include <memory>
#include <limits>
#include <atomic>
#include "seal/util/globals.h"
#include "seal/util/common.h"
#include "seal/util/locks.h"

namespace seal
{
    namespace util
    {
        // Largest size of single allocation that can be requested from memory pool
        constexpr std::int64_t mempool_max_single_alloc_uint64_count = (1LL << 48) - 1;

        // Number of different size allocations allowed by a single memory pool
        constexpr std::int64_t mempool_max_pool_head_count = 
            std::numeric_limits<std::int64_t>::max();

        // Largest allowed size of batch allocation
        constexpr std::int64_t mempool_max_batch_alloc_uint64_count = (1LL << 48) - 1;

        constexpr std::int64_t mempool_first_alloc_count = 1;

        constexpr double mempool_alloc_size_multiplier = 1.05;

        class MemoryPoolItem
        {
        public:
            MemoryPoolItem(std::uint64_t *data) : data_(data), next_(nullptr)
            {
            }

            inline std::uint64_t *data()
            {
                return data_;
            }

            inline const std::uint64_t *data() const
            {
                return data_;
            }

            inline MemoryPoolItem* &next()
            {
                return next_;
            }

            inline const MemoryPoolItem *next() const
            {
                return next_;
            }

        private:
            MemoryPoolItem(const MemoryPoolItem &copy) = delete;

            MemoryPoolItem &operator =(const MemoryPoolItem &assign) = delete;

            std::uint64_t *data_;

            MemoryPoolItem *next_;
        };

        class MemoryPoolHead
        {
        public:
            struct allocation
            {
                allocation() : 
                    size(0), data_ptr(nullptr), free(0), head_ptr(nullptr)
                {
                }

                // Size of an allocation (in uint64_count)
                std::int64_t size;

                // Pointer to start of allocation
                std::uint64_t *data_ptr;

                // How much free space is left (in uint64_count)
                std::int64_t free;

                // Pointer to current head of allocation
                std::uint64_t *head_ptr;
            };

            virtual ~MemoryPoolHead()
            {
            }

            virtual std::int64_t uint64_count() const = 0;

            virtual std::int64_t alloc_item_count() const = 0;

            virtual MemoryPoolItem *get() = 0;

            virtual void add(MemoryPoolItem *new_first) = 0;
        };

        class MemoryPoolHeadMT : public MemoryPoolHead
        {
        public:
            // Creates a new MemoryPoolHeadMT with allocation for one single item.
            MemoryPoolHeadMT(std::int64_t uint64_count);

            ~MemoryPoolHeadMT() override;

            inline std::int64_t uint64_count() const override
            {
                return uint64_count_;
            }

            // Returns the total number of items allocated
            inline std::int64_t alloc_item_count() const override
            {
                return alloc_item_count_;
            }

            MemoryPoolItem *get() override;

            inline void add(MemoryPoolItem *new_first) override
            {
                bool expected = false;
                while (!locked_.compare_exchange_strong(expected, true, std::memory_order_acquire))
                {
                    expected = false;
                }
                MemoryPoolItem *old_first = first_item_;
                new_first->next() = old_first;
                first_item_ = new_first;
                locked_.store(false, std::memory_order_release);
            }

        private:
            MemoryPoolHeadMT(const MemoryPoolHeadMT &copy) = delete;

            MemoryPoolHeadMT &operator =(const MemoryPoolHeadMT &assign) = delete;

            mutable std::atomic<bool> locked_;

            const std::int64_t uint64_count_;

            volatile std::int64_t alloc_item_count_;

            std::vector<allocation> allocs_;

            MemoryPoolItem* volatile first_item_;
        };

        class MemoryPoolHeadST : public MemoryPoolHead
        {
        public:
            // Creates a new MemoryPoolHeadST with allocation for one single item.
            MemoryPoolHeadST(std::int64_t uint64_count);

            ~MemoryPoolHeadST() override;

            inline std::int64_t uint64_count() const override
            {
                return uint64_count_;
            }

            // Returns the total number of items allocated
            inline std::int64_t alloc_item_count() const override
            {
                return alloc_item_count_;
            }

            MemoryPoolItem *get() override;

            inline void add(MemoryPoolItem *new_first) override
            {
                new_first->next() = first_item_;
                first_item_ = new_first;
            }

        private:
            MemoryPoolHeadST(const MemoryPoolHeadST &copy) = delete;

            MemoryPoolHeadST &operator =(const MemoryPoolHeadST &assign) = delete;

            std::int64_t uint64_count_;

            std::int64_t alloc_item_count_;

            std::vector<allocation> allocs_;

            MemoryPoolItem *first_item_;
        };

        class ConstPointer;

        class Pointer
        {
        public:
            friend class ConstPointer;

            Pointer() : data_(nullptr), head_(nullptr), item_(nullptr), alias_(false)
            {
            }

            Pointer(MemoryPoolHead *head) : data_(nullptr), head_(nullptr), item_(nullptr), 
                alias_(false)
            {
#ifdef SEAL_DEBUG
                if (head == nullptr)
                {
                    throw std::invalid_argument("head");
                }
#endif
                head_ = head;
                item_ = head->get();
                data_ = item_->data();
            }

            Pointer(Pointer &&move) noexcept : data_(move.data_), head_(move.head_), 
                item_(move.item_), alias_(move.alias_)
            {
                move.data_ = nullptr;
                move.head_ = nullptr;
                move.item_ = nullptr;
                move.alias_ = false;
            }

            inline std::uint64_t &operator [](std::size_t index)
            {
                return data_[index];
            }

            inline const std::uint64_t &operator [](std::size_t index) const
            {
                return data_[index];
            }

            inline Pointer &operator =(Pointer &&assign)
            {
                acquire(assign);
                return *this;
            }

            inline bool is_set() const
            {
                return data_ != nullptr;
            }

            inline std::uint64_t *get()
            {
                return data_;
            }

            inline const std::uint64_t *get() const
            {
                return data_;
            }

            inline bool is_alias() const
            {
                return alias_;
            }

            inline void release()
            {
                if (head_ != nullptr)
                {
                    head_->add(item_);
                }
                else if (data_ != nullptr && !alias_)
                {
                    delete[] data_;
                }
                data_ = nullptr;
                head_ = nullptr;
                item_ = nullptr;
                alias_ = false;
            }

            inline void acquire(Pointer &other)
            {
                if (this == &other)
                {
                    return;
                }
                if (head_ != nullptr)
                {
                    head_->add(item_);
                }
                else if (data_ != nullptr && !alias_)
                {
                    delete[] data_;
                }
                data_ = other.data_;
                head_ = other.head_;
                item_ = other.item_;
                alias_ = other.alias_;
                other.data_ = nullptr;
                other.head_ = nullptr;
                other.item_ = nullptr;
                other.alias_ = false;
            }

            inline void swap_with(Pointer &other) noexcept 
            {
                std::swap(data_, other.data_);
                std::swap(head_, other.head_);
                std::swap(item_, other.item_);
                std::swap(alias_, other.alias_);
            }

            inline void swap_with(Pointer &&other) noexcept
            {
                std::swap(data_, other.data_);
                std::swap(head_, other.head_);
                std::swap(item_, other.item_);
                std::swap(alias_, other.alias_);
            }
            
            ~Pointer()
            {
                release();
            }

            inline static Pointer Owning(std::uint64_t *pointer)
            {
                return Pointer(pointer, false);
            }

            inline static Pointer Aliasing(std::uint64_t *pointer)
            {
                return Pointer(pointer, true);
            }

            Pointer(std::uint64_t *pointer, bool alias) noexcept : data_(pointer), head_(nullptr), 
                item_(nullptr), alias_(alias)
            {
            }

        private:
            Pointer(const Pointer &copy) = delete;

            Pointer &operator =(const Pointer &assign) = delete;

            std::uint64_t *data_;

            MemoryPoolHead *head_;

            MemoryPoolItem *item_;

            bool alias_;
        };

        class ConstPointer
        {
        public:
            ConstPointer() : data_(nullptr), head_(nullptr), item_(nullptr), alias_(false)
            {
            }

            ConstPointer(MemoryPoolHead *head) : data_(nullptr), head_(nullptr), item_(nullptr), 
                alias_(false)
            {
#ifdef SEAL_DEBUG
                if (head == nullptr)
                {
                    throw std::invalid_argument("head");
                }
#endif
                data_ = item_->data();
                head_ = head;
                item_ = head->get();
            }

            ConstPointer(ConstPointer &&move) noexcept : data_(move.data_), head_(move.head_), 
                item_(move.item_), alias_(move.alias_)
            {
                move.data_ = nullptr;
                move.head_ = nullptr;
                move.item_ = nullptr;
                move.alias_ = false;
            }

            ConstPointer(Pointer &&move) noexcept : data_(move.data_), head_(move.head_), 
                item_(move.item_), alias_(move.alias_)
            {
                move.data_ = nullptr;
                move.head_ = nullptr;
                move.item_ = nullptr;
                move.alias_ = false;
            }

            inline ConstPointer &operator =(ConstPointer &&assign)
            {
                acquire(assign);
                return *this;
            }

            inline ConstPointer &operator =(Pointer &&assign)
            {
                acquire(assign);
                return *this;
            }

            inline std::uint64_t operator [](int index) const
            {
                return data_[index];
            }

            inline bool is_set() const
            {
                return data_ != nullptr;
            }

            inline const std::uint64_t *get() const
            {
                return data_;
            }

            inline void release()
            {
                if (head_ != nullptr)
                {
                    head_->add(item_);
                }
                else if (data_ != nullptr && !alias_)
                {
                    delete[] data_;
                }
                data_ = nullptr;
                head_ = nullptr;
                item_ = nullptr;
                alias_ = false;
            }

            inline void acquire(ConstPointer &other)
            {
                if (this == &other)
                {
                    return;
                }
                if (head_ != nullptr)
                {
                    head_->add(item_);
                }
                else if (data_ != nullptr && !alias_)
                {
                    delete[] data_;
                }                
                data_ = other.data_;
                head_ = other.head_;
                item_ = other.item_;
                alias_ = other.alias_;
                other.data_ = nullptr;
                other.head_ = nullptr;
                other.item_ = nullptr;
                other.alias_ = false;
            }

            inline void acquire(Pointer &other)
            {
                if (head_ != nullptr)
                {
                    head_->add(item_);
                }
                else if (data_ != nullptr && !alias_)
                {
                    delete[] data_;
                }                
                data_ = other.data_;
                head_ = other.head_;
                item_ = other.item_;
                alias_ = other.alias_;
                other.data_ = nullptr;
                other.head_ = nullptr;
                other.item_ = nullptr;
                other.alias_ = false;
            }

            inline void swap_with(ConstPointer &other) noexcept
            {
                std::swap(data_, other.data_);
                std::swap(head_, other.head_);
                std::swap(item_, other.item_);
                std::swap(alias_, other.alias_);
            }

            inline void swap_with(ConstPointer &&other) noexcept
            {
                std::swap(data_, other.data_);
                std::swap(head_, other.head_);
                std::swap(item_, other.item_);
                std::swap(alias_, other.alias_);
            }

            ~ConstPointer()
            {
                release();
            }

            inline static ConstPointer Owning(std::uint64_t *pointer)
            {
                return ConstPointer(pointer, false);
            }

            inline static ConstPointer Aliasing(const std::uint64_t *pointer)
            {
                return ConstPointer(const_cast<uint64_t*>(pointer), true);
            }

        private:
            ConstPointer(std::uint64_t *pointer, bool alias) : data_(pointer), head_(nullptr), 
                item_(nullptr), alias_(alias)
            {
            }

            ConstPointer(const ConstPointer &copy) = delete;

            ConstPointer &operator =(const ConstPointer &assign) = delete;

            std::uint64_t *data_;

            MemoryPoolHead *head_;

            MemoryPoolItem *item_;

            bool alias_;
        };

        class MemoryPool
        {
        public:
            virtual ~MemoryPool()
            {
            }

            virtual Pointer get_for_uint64_count(std::int64_t uint64_count) = 0;

            virtual std::int64_t pool_count() const = 0;

            virtual std::int64_t alloc_uint64_count() const = 0;

            virtual std::int64_t alloc_byte_count() const = 0;
        };

        class MemoryPoolMT : public MemoryPool
        {
        public:
            MemoryPoolMT() 
            {
            }

            ~MemoryPoolMT() override;

            Pointer get_for_uint64_count(std::int64_t uint64_count) override;

            inline std::int64_t pool_count() const override
            {
                ReaderLock lock(pools_locker_.acquire_read());
                return static_cast<std::int64_t>(pools_.size());
            }

            std::int64_t alloc_uint64_count() const override;

            std::int64_t alloc_byte_count() const override
            {
                return alloc_uint64_count() * bytes_per_uint64;
            }

        private:
            MemoryPoolMT(const MemoryPoolMT &copy) = delete;

            MemoryPoolMT &operator =(const MemoryPoolMT &assign) = delete;

            mutable ReaderWriterLocker pools_locker_;

            std::vector<MemoryPoolHead*> pools_;
        };

        class MemoryPoolST : public MemoryPool
        {
        public:
            MemoryPoolST()
            {
            }

            ~MemoryPoolST() override;

            Pointer get_for_uint64_count(std::int64_t uint64_count) override;

            inline std::int64_t pool_count() const override
            {
                return static_cast<std::int64_t>(pools_.size());
            }

            std::int64_t alloc_uint64_count() const override;
            
            std::int64_t alloc_byte_count() const override
            {
                return alloc_uint64_count() * bytes_per_uint64;
            }

        private:
            MemoryPoolST(const MemoryPoolST &copy) = delete;

            MemoryPoolST &operator =(const MemoryPoolST &assign) = delete;

            std::vector<MemoryPoolHead*> pools_;
        };

        inline Pointer duplicate_if_needed(std::uint64_t *original, 
            std::int64_t uint64_count, bool condition, MemoryPool &pool)
        {
#ifdef SEAL_DEBUG
            if (original == nullptr && uint64_count > 0)
            {
                throw std::invalid_argument("original");
            }
            if (uint64_count < 0)
            {
                throw std::invalid_argument("uint64_count");
            }
#endif
            if (condition == false)
            {
                return Pointer::Aliasing(original);
            }
            if (!uint64_count)
            {
                return Pointer();
            }
            Pointer allocation(pool.get_for_uint64_count(uint64_count));
            std::memcpy(allocation.get(), original, 
                static_cast<std::size_t>(uint64_count) * bytes_per_uint64);
            return allocation;
        }

        inline ConstPointer duplicate_if_needed(const std::uint64_t *original, 
            std::int64_t uint64_count, bool condition, MemoryPool &pool)
        {
#ifdef SEAL_DEBUG
            if (original == nullptr && uint64_count > 0)
            {
                throw std::invalid_argument("original");
            }
            if (uint64_count < 0)
            {
                throw std::invalid_argument("uint64_count");
            }
#endif
            if (condition == false)
            {
                return ConstPointer::Aliasing(original);
            }
            if (!uint64_count)
            {
                return ConstPointer();
            }
            Pointer allocation = pool.get_for_uint64_count(uint64_count);
            std::memcpy(allocation.get(), original, 
                static_cast<std::size_t>(uint64_count) * bytes_per_uint64);
            ConstPointer const_allocation;
            const_allocation.acquire(allocation);
            return const_allocation;
        }
    }
}
