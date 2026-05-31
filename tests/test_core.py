import unittest

from pylinearcode import (
    GF,
    bch_code,
    bch_generator_polynomial,
    extended_hamming_code,
    hamming_code,
    macwilliams_transform,
    parity_check_code,
    reed_muller_code,
    reed_solomon_code,
    repetition_code,
    simplex_code,
)
from pylinearcode.polynomials import divmod, evaluate, gcd, mul


class FiniteFieldTests(unittest.TestCase):
    def test_gf4_arithmetic(self):
        field = GF(4)
        alpha = 2
        self.assertEqual(field.mul(alpha, alpha), field.add(alpha, 1))
        self.assertEqual(field.pow(alpha, 3), 1)
        self.assertEqual(field.div(alpha, alpha), 1)


class PolynomialTests(unittest.TestCase):
    def test_polynomial_arithmetic_over_gf2(self):
        field = GF(2)
        product = mul(field, [1, 1], [1, 0, 1])
        self.assertEqual(product, [1, 1, 1, 1])
        quotient, remainder = divmod(field, product, [1, 1])
        self.assertEqual(quotient, [1, 0, 1])
        self.assertEqual(remainder, [0])
        self.assertEqual(gcd(field, product, [1, 1]), [1, 1])
        self.assertEqual(evaluate(field, [1, 1, 1], 1), 1)


class LinearCodeTests(unittest.TestCase):
    def test_binary_hamming_code_parameters_and_decode(self):
        field = GF(2)
        code = hamming_code(field, 3)
        self.assertEqual((code.length, code.dimension), (7, 4))
        self.assertEqual(code.minimum_distance(), 3)

        word = code.encode([1, 0, 1, 1])
        received = word[:]
        received[2] ^= 1
        decoded, error = code.syndrome_decode(received, max_weight=1)
        self.assertEqual(decoded, word)
        self.assertEqual(error[2], 1)

    def test_extended_hamming_and_simplex(self):
        field = GF(2)
        extended = extended_hamming_code(field, 3)
        simplex = simplex_code(field, 3)
        self.assertEqual((extended.length, extended.dimension), (8, 4))
        self.assertEqual(extended.minimum_distance(), 4)
        self.assertEqual(simplex.weight_distribution()[4], 7)

    def test_reed_solomon_is_mds(self):
        field = GF(5)
        code = reed_solomon_code(field, length=5, dimension=3)
        self.assertEqual(code.minimum_distance(), 3)
        self.assertTrue(code.is_mds())

    def test_repetition_dual_is_parity_check(self):
        field = GF(2)
        repetition = repetition_code(field, 5)
        parity = parity_check_code(field, 5)
        self.assertEqual(repetition.dual().dimension, parity.dimension)
        for row in repetition.dual().generator_matrix:
            self.assertTrue(parity.contains(row))

    def test_shorten_and_puncture(self):
        field = GF(2)
        code = hamming_code(field, 3)
        punctured = code.puncture([0])
        shortened = code.shorten([0])
        self.assertEqual(punctured.length, 6)
        self.assertEqual(shortened.length, 6)
        self.assertEqual(shortened.dimension, 3)

    def test_hull(self):
        field = GF(2)
        code = hamming_code(field, 3)
        hull = code.hull()
        self.assertTrue(all(code.contains(row) for row in hull.generator_matrix))
        self.assertTrue(all(code.dual().contains(row) for row in hull.generator_matrix))

    def test_bch_code(self):
        field = GF(2)
        self.assertEqual(bch_generator_polynomial(field, extension_degree=3, designed_distance=3), [1, 0, 1, 1])
        code = bch_code(field, extension_degree=4, designed_distance=5)
        self.assertEqual((code.length, code.dimension), (15, 7))
        self.assertEqual(code.minimum_distance(), 5)

    def test_reed_muller_code(self):
        code = reed_muller_code(1, 3)
        self.assertEqual((code.length, code.dimension), (8, 4))
        self.assertEqual(code.minimum_distance(), 4)
        self.assertEqual(code.weight_distribution()[4], 14)

    def test_dual_weight_distribution_and_covering_radius(self):
        field = GF(2)
        code = hamming_code(field, 3)
        dual_distribution = code.dual_weight_distribution()
        self.assertEqual(dual_distribution[0], 1)
        self.assertEqual(dual_distribution[4], 7)
        self.assertEqual(code.covering_radius(max_weight=1), 1)
        self.assertEqual(
            macwilliams_transform(code.weight_distribution(), q=2, length=7, dimension=4),
            dual_distribution,
        )


if __name__ == "__main__":
    unittest.main()
