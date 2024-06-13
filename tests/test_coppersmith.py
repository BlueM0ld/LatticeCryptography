import pytest
from src import coppersmith_temp
from sage.all import Integer, random_prime


def test_coppersmith_univariate():
    # RSA modulus N
    N = random_prime(2**150) * random_prime(2**150)
    message = Integer('thepasswordfortodayisswordfish', base=35)

    # Encrypt with e=3
    e = 3
    c = message**e % N

    # Define the known and unknown parts of the message
    known_prefix = 'thepasswordfortodayis000000000'
    unknown_length = 9  # Length of the unknown part

    Q, recovered_root, recovered_message = coppersmith_temp.coppersmith_univariate(
        N, c, known_prefix, unknown_length, e)

    assert recovered_root is not None, "Failed to recover the root"
    assert recovered_message == 'swordfish', "Recovered message is incorrect"


if __name__ == '__main__':
    pytest.main()
