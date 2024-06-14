import pytest
from src import coppersmith_univariate
from sage.all import Integer, random_prime


def setup_rsa_and_cipher(known_prefix, unknown_length, e):
    N = random_prime(2**150) * random_prime(2**150)
    message = Integer('thepasswordfortodayisswordfish', base=35)
    c = message**e % N
    return N, c, known_prefix, unknown_length, e


def test_coppersmith_univariate_case1():
    N, c, known_prefix, unknown_length, e = setup_rsa_and_cipher(
        'thepasswordfortodayis000000000', 9, 3)
    recovered_message = coppersmith_univariate.coppersmith_univariate(
        N, c, known_prefix, unknown_length, e)
    print(recovered_message)


def test_coppersmith_univariate_case2():
    N, c, known_prefix, unknown_length, e = setup_rsa_and_cipher(
        'thepasswordfortodayiss00000000', 8, 3)
    recovered_message = coppersmith_univariate.coppersmith_univariate(
        N, c, known_prefix, unknown_length, e)
    print(recovered_message)


def test_coppersmith_univariate_case3():
    N, c, known_prefix, unknown_length, e = setup_rsa_and_cipher(
        'thepasswordfortodayis000000000', 9, 5)
    recovered_message = coppersmith_univariate.coppersmith_univariate(
        N, c, known_prefix, unknown_length, e)
    print(recovered_message)


def test_coppersmith_univariate_case4():
    N, c, known_prefix, unknown_length, e = setup_rsa_and_cipher(
        'thepasswordfortodayiss00000000', 8, 5)
    recovered_message = coppersmith_univariate.coppersmith_univariate(
        N, c, known_prefix, unknown_length, e)
    print(recovered_message)


def test_coppersmith_univariate_case5():
    N, c, known_prefix, unknown_length, e = setup_rsa_and_cipher(
        'thepasswordfortodayis000000000', 9, 7)
    recovered_message = coppersmith_univariate.coppersmith_univariate(
        N, c, known_prefix, unknown_length, e)
    print(recovered_message)


def test_coppersmith_univariate_case6():
    N, c, known_prefix, unknown_length, e = setup_rsa_and_cipher(
        'thepasswordfortodayiss00000000', 8, 7)
    recovered_message = coppersmith_univariate.coppersmith_univariate(
        N, c, known_prefix, unknown_length, e)
    print(recovered_message)


def test_coppersmith_univariate_case7():
    N, c, known_prefix, unknown_length, e = setup_rsa_and_cipher(
        'thepasswordfortodayis000000000', 9, 11)
    recovered_message = coppersmith_univariate.coppersmith_univariate(
        N, c, known_prefix, unknown_length, e)
    print(recovered_message)


def test_coppersmith_univariate_case8():
    N, c, known_prefix, unknown_length, e = setup_rsa_and_cipher(
        'thepasswordfortodayiss00000000', 8, 13)
    recovered_message = coppersmith_univariate.coppersmith_univariate(
        N, c, known_prefix, unknown_length, e)
    print(recovered_message)


def test_coppersmith_univariate_case9():
    N, c, known_prefix, unknown_length, e = setup_rsa_and_cipher(
        'thepasswordfortodayis000000000', 9, 19)
    recovered_message = coppersmith_univariate.coppersmith_univariate(
        N, c, known_prefix, unknown_length, e)
    print(recovered_message)


def test_coppersmith_univariate_case10():
    N, c, known_prefix, unknown_length, e = setup_rsa_and_cipher(
        'thepasswordfortodayiss00000000', 8, 23)
    recovered_message = coppersmith_univariate.coppersmith_univariate(
        N, c, known_prefix, unknown_length, e)
    print(recovered_message)


if __name__ == '__main__':
    pytest.main()
