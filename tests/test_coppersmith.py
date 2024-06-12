import unittest
from sage.all import Integer, random_prime
from coppersmith_temp import coppersmith_univariate


class TestCoppersmithUnivariate(unittest.TestCase):
    def test_coppersmith_univariate(self):
        # RSA modulus N
        N = random_prime(2**150) * random_prime(2**150)
        message = Integer('thepasswordfortodayisswordfish', base=35)

        # Encrypt with e=3
        e = 3
        c = message**e % N

        # Define the known and unknown parts of the message
        known_prefix = 'thepasswordfortodayis000000000'
        unknown_length = 9  # Length of the unknown part

        Q, recovered_root, recovered_message = coppersmith_univariate(
            N, c, known_prefix, unknown_length, e)

        self.assertIsNotNone(recovered_root, "Failed to recover the root")
        self.assertEqual(recovered_message, 'swordfish',
                         "Recovered message is incorrect")


if __name__ == '__main__':
    unittest.main()
