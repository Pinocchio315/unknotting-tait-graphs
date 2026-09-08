"""Check the prime-by-prime bound independently of any deposited pair count."""
import unittest
from connected_sum_generators import prime_ranks, sharp_primes


class ConnectedSumGeneratorTests(unittest.TestCase):
    def test_coprime_cyclic_summands_do_not_add_generator_numbers(self):
        ranks = prime_ranks([3, 5])
        self.assertEqual(ranks, {3: 1, 5: 1})
        self.assertEqual(max(ranks.values()), 1)
        self.assertFalse(sharp_primes([3], 1) & sharp_primes([5], 1))

    def test_sharp_u_two_bounds_require_a_common_prime(self):
        # Double-cover groups of 8_18 and 9_40 respectively. Each needs two
        # generators, but their direct sum needs three, rather than four.
        self.assertEqual(sharp_primes([3, 15], 2), {3})
        self.assertEqual(sharp_primes([5, 15], 2), {5})
        self.assertFalse(sharp_primes([3, 15], 2) & sharp_primes([5, 15], 2))
        self.assertEqual(prime_ranks([3, 15, 5, 15]), {3: 3, 5: 3})

    def test_common_prime_does_certify_addition(self):
        self.assertEqual(sharp_primes([3, 15], 2) & sharp_primes([3, 21], 2), {3})
        self.assertEqual(prime_ranks([3, 15, 3, 21])[3], 4)
        # Repeated powers in one cyclic factor do not count as extra generators.
        self.assertEqual(prime_ranks([9, 27, 25]), {3: 2, 5: 1})

    def test_invalid_or_contradictory_exact_inputs_are_rejected(self):
        for factors in ([0], [-3], [True]):
            with self.subTest(factors=factors), self.assertRaises(ValueError):
                prime_ranks(factors)
        with self.assertRaises(ValueError):
            sharp_primes([3, 3], 1)


if __name__ == '__main__':
    unittest.main()
