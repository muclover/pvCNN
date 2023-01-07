#pragma once

#include "seal/util/defines.h"

#ifdef SEAL_USE_SHARED_MUTEX_FOR_RW_LOCK
#include <shared_mutex>
#include <mutex>

namespace seal
{
    namespace util
    {
        using ReaderLock = std::shared_lock<std::shared_mutex>;

        using WriterLock = std::unique_lock<std::shared_mutex>;

        class ReaderWriterLocker
        {
        public:
            ReaderWriterLocker() = default;

            inline ReaderLock acquire_read()
            {
                return ReaderLock(rw_lock_mutex_);
            }

            inline WriterLock acquire_write()
            {
                return WriterLock(rw_lock_mutex_);
            }

            inline ReaderLock try_acquire_read()
            {
                return ReaderLock(rw_lock_mutex_, std::try_to_lock);
            }

            inline WriterLock try_acquire_write()
            {
                return WriterLock(rw_lock_mutex_, std::try_to_lock);
            }

        private:
            ReaderWriterLocker(const ReaderWriterLocker &copy) = delete;

            ReaderWriterLocker &operator =(const ReaderWriterLocker &assign) = delete;

            std::shared_mutex rw_lock_mutex_;
        };
    }
}
#else
#include <atomic>

namespace seal
{
    namespace util
    {
        struct try_to_lock_t
        {
        };

        constexpr try_to_lock_t try_to_lock{};

        class ReaderWriterLocker;

        class ReaderLock
        {
        public:
            ReaderLock() : locker_(nullptr)
            {
            }

            ReaderLock(ReaderLock &&move) noexcept : locker_(move.locker_)
            {
                move.locker_ = nullptr;
            }

            ReaderLock(ReaderWriterLocker &locker) : locker_(nullptr)
            {
                acquire(locker);
            }

            ReaderLock(ReaderWriterLocker &locker, try_to_lock_t) : locker_(nullptr)
            {
                try_acquire(locker);
            }

            ~ReaderLock()
            {
                unlock();
            }

            inline bool owns_lock() const
            {
                return locker_ != nullptr;
            }

            void unlock();

            inline void swap_with(ReaderLock &lock)
            {
                std::swap(locker_, lock.locker_);
            }

            inline ReaderLock &operator =(ReaderLock &&lock)
            {
                swap_with(lock);
                lock.unlock();
                return *this;
            }

        private:
            void acquire(ReaderWriterLocker &locker);

            bool try_acquire(ReaderWriterLocker &locker);

            ReaderWriterLocker *locker_;
        };

        class WriterLock
        {
        public:
            WriterLock() : locker_(nullptr)
            {
            }

            WriterLock(WriterLock &&move) noexcept : locker_(move.locker_)
            {
                move.locker_ = nullptr;
            }

            WriterLock(ReaderWriterLocker &locker) : locker_(nullptr)
            {
                acquire(locker);
            }

            WriterLock(ReaderWriterLocker &locker, try_to_lock_t) : locker_(nullptr)
            {
                try_acquire(locker);
            }

            ~WriterLock()
            {
                unlock();
            }

            inline bool owns_lock() const
            {
                return locker_ != nullptr;
            }

            void unlock();

            inline void swap_with(WriterLock &lock)
            {
                std::swap(locker_, lock.locker_);
            }

            inline WriterLock &operator =(WriterLock &&lock)
            {
                swap_with(lock);
                lock.unlock();
                return *this;
            }

        private:
            void acquire(ReaderWriterLocker &locker);

            bool try_acquire(ReaderWriterLocker &locker);

            ReaderWriterLocker *locker_;
        };

        class ReaderWriterLocker
        {
            friend class ReaderLock;

            friend class WriterLock;

        public:
            ReaderWriterLocker() : reader_locks_(0), writer_locked_(false)
            {
            }

            inline ReaderLock acquire_read()
            {
                return ReaderLock(*this);
            }

            inline WriterLock acquire_write()
            {
                return WriterLock(*this);
            }

            inline ReaderLock try_acquire_read()
            {
                return ReaderLock(*this, try_to_lock);
            }

            inline WriterLock try_acquire_write()
            {
                return WriterLock(*this, try_to_lock);
            }

        private:
            ReaderWriterLocker(const ReaderWriterLocker &copy) = delete;

            ReaderWriterLocker &operator =(const ReaderWriterLocker &assign) = delete;

            std::atomic<int> reader_locks_;

            std::atomic<bool> writer_locked_;
        };

        inline void ReaderLock::unlock()
        {
            if (locker_ == nullptr)
            {
                return;
            }
            locker_->reader_locks_.fetch_sub(1, std::memory_order_release);
            locker_ = nullptr;
        }

        inline void ReaderLock::acquire(ReaderWriterLocker &locker)
        {
            unlock();
            do
            {
                locker.reader_locks_.fetch_add(1, std::memory_order_acquire);
                locker_ = &locker;
                if (locker.writer_locked_.load(std::memory_order_acquire))
                {
                    unlock();
                    while (locker.writer_locked_.load(std::memory_order_acquire));
                }
            } while (locker_ == nullptr);
        }

        inline bool ReaderLock::try_acquire(ReaderWriterLocker &locker)
        {
            unlock();
            locker.reader_locks_.fetch_add(1, std::memory_order_acquire);
            locker_ = &locker;
            if (locker.writer_locked_.load(std::memory_order_acquire))
            {
                unlock();
                return false;
            }
            return true;
        }

        inline void WriterLock::acquire(ReaderWriterLocker &locker)
        {
            unlock();
            bool expected = false;
            while (!locker.writer_locked_.compare_exchange_strong(expected, true, std::memory_order_acquire))
            {
                expected = false;
            }
            locker_ = &locker;
            while (locker.reader_locks_.load(std::memory_order_acquire) != 0);
        }

        inline bool WriterLock::try_acquire(ReaderWriterLocker &locker)
        {
            unlock();
            bool expected = false;
            if (!locker.writer_locked_.compare_exchange_strong(expected, true, std::memory_order_acquire))
            {
                return false;
            }
            locker_ = &locker;
            if (locker.reader_locks_.load(std::memory_order_acquire) != 0)
            {
                unlock();
                return false;
            }
            return true;
        }

        inline void WriterLock::unlock()
        {
            if (locker_ == nullptr)
            {
                return;
            }
            locker_->writer_locked_.store(false, std::memory_order_release);
            locker_ = nullptr;
        }
    }
}
#endif
