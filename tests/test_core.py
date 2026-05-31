import unittest

from pylinearcode import GF, hamming_code, parity_check_code, reed_solomon_code, repetition_code


class FiniteFieldTests(unittest.TestCase):
    def test_gf4_arithmetic(self):
        field = GF(4)
        alpha = 2
        self.assertEqual(field.mul(alpha, alpha), field.add(alpha, 1))
        self.assertEqual(field.pow(alpha, 3), 1)
        self.assertEqual(field.div(alpha, alpha), 1)


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


if __name__ == "__main__":
    unittest.main()

