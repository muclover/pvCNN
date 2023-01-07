#pragma once

#include <memory>
#include <utility>
#include "seal/context.h"
#include "seal/util/polymodulus.h"
#include "seal/util/smallntt.h"
#include "seal/memorypoolhandle.h"
#include "seal/publickey.h"
#include "seal/secretkey.h"
#include "seal/evaluationkeys.h"
#include "seal/galoiskeys.h"

namespace seal
{
    /**
    Generates matching secret key and public key. An existing KeyGenerator can also at any time
    be used to generate evaluation keys and Galois keys. Constructing a KeyGenerator requires 
    only a SEALContext.

    @see EncryptionParameters for more details on encryption parameters.
    @see SecretKey for more details on secret key.
    @see PublicKey for more details on public key.
    @see EvaluationKeys for more details on evaluation keys.
    @see GaloisKeys for more details on Galois keys.
    */
    class KeyGenerator
    {
    public:
        /**
        Creates a KeyGenerator initialized with the specified SEALContext.

        @param[in] context The SEALContext
        @throws std::invalid_argument if encryption parameters is not valid
        */
        KeyGenerator(const SEALContext &context);

        /**
        Creates an KeyGenerator instance initialized with the specified SEALContext and 
        specified previously secret and public keys. This can e.g. be used to increase the 
        number of evaluation keys from what had earlier been generated, or to generate 
        Galois keys in case they had not been generated earlier.

        @param[in] context The SEALContext
        @param[in] secret_key A previously generated secret key
        @param[in] public_key A previously generated public key
        @throws std::invalid_argument if encryption parameters are not valid
        @throws std::invalid_argument if secret_key or public_key is not valid for 
        encryption parameters
        */
        KeyGenerator(const SEALContext &context, const SecretKey &secret_key, 
            const PublicKey &public_key);

        /**
        Returns a const reference to the secret key.
        */
        const SecretKey &secret_key() const;

        /**
        Returns a const reference to the public key.
        */
        const PublicKey &public_key() const;

        /**
        Generates the specified number of evaluation keys.

        @param[in] decomposition_bit_count The decomposition bit count
        @param[in] count The number of evaluation keys to generate
        @param[out] evaluation_keys The evaluation keys instance to overwrite with the 
        generated keys
        @throws std::invalid_argument if decomposition_bit_count is not within [1, 60]
        @throws std::invalid_argument if count is negative
        */
        void generate_evaluation_keys(int decomposition_bit_count, int count, 
            EvaluationKeys &evaluation_keys);

        /**
        Generates evaluation keys containing one key.

        @param[in] decomposition_bit_count The decomposition bit count
        @param[out] evaluation_keys The evaluation keys instance to overwrite with the
        generated keys
        @throws std::invalid_argument if decomposition_bit_count is not within [1, 60]
        */
        inline void generate_evaluation_keys(int decomposition_bit_count, 
            EvaluationKeys &evaluation_keys)
        {
            generate_evaluation_keys(decomposition_bit_count, 1, evaluation_keys);
        }

        /**
        Generates Galois keys. This function creates logarithmically many (in degree of the
        polynomial modulus) Galois keys that is sufficient to apply any Galois automorphism
        (e.g. rotations) on encrypted data. Most users will want to use this overload of
        the function.

        @param[in] decomposition_bit_count The decomposition bit count
        @param[out] galois_keys The Galois keys instance to overwrite with the generated keys
        @throws std::invalid_argument if decomposition_bit_count is not within [1, 60]
        */
        void generate_galois_keys(int decomposition_bit_count, GaloisKeys &galois_keys);

        /**
        Generates Galois keys. This function creates specific Galois keys that can be used to
        apply specific Galois automorphisms on encrypted data. The user needs to give as 
        input a vector of Galois elements corresponding to the keys that are to be created.
        
        The Galois elements are odd integers in the interval [1, M-1], where M = 2*N, and
        N = degree(poly_modulus). Used with batching, a Galois element 3^i % M corresponds
        to a cyclic row rotation i steps to the left, and a Galois element 3^(N/2-i) % M
        corresponds to a cyclic row rotation i steps to the right. The Galois element M-1
        corresponds to a column rotation (row swap). In the polynomial view (not batching),
        a Galois automorphism by a Galois element p changes Enc(plain(x)) to Enc(plain(x^p)).

        @param[in] decomposition_bit_count The decomposition bit count
        @param[in] galois_elts The Galois elements for which to generate keys
        @param[out] galois_keys The Galois keys instance to overwrite with the generated keys
        @throws std::invalid_argument if decomposition_bit_count is not within [1, 60]
        @throws std::invalid_argument if the Galois elements are not valid
        */
        void generate_galois_keys(int decomposition_bit_count,
            const std::vector<std::uint64_t> &galois_elts, GaloisKeys &galois_keys);

    private:
        KeyGenerator(const KeyGenerator &copy) = delete;

        KeyGenerator &operator =(const KeyGenerator &assign) = delete;

        KeyGenerator(KeyGenerator &&source) = delete;

        KeyGenerator &operator =(KeyGenerator &&assign) = delete;

        void set_poly_coeffs_zero_one_negone(std::uint64_t *poly, UniformRandomGenerator *random) const;

        void set_poly_coeffs_normal(std::uint64_t *poly, UniformRandomGenerator *random) const;

        void set_poly_coeffs_uniform(std::uint64_t *poly, UniformRandomGenerator *random);

        void compute_secret_key_array(int max_power);

        void populate_decomposition_factors(int decomposition_bit_count, 
            std::vector<std::vector<std::uint64_t> > &decomposition_factors);

        /**
        Generates new matching set of secret key and public key.
        */
        void generate();

        /**
        Returns whether secret key and public key have been generated.
        */
        inline bool is_generated() const
        {
            return generated_;
        }

        MemoryPoolHandle pool_ = MemoryPoolHandle::Global();

        EncryptionParameters parms_;

        EncryptionParameterQualifiers qualifiers_;

        std::vector<util::SmallNTTTables> small_ntt_tables_;

        PublicKey public_key_;

        SecretKey secret_key_;

        UniformRandomGeneratorFactory *random_generator_ = nullptr;

        util::PolyModulus polymod_;

        int secret_key_array_size_;

        util::Pointer secret_key_array_;

        mutable util::ReaderWriterLocker secret_key_array_locker_;

        bool generated_ = false;
    };
}
