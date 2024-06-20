import pytest
from src import coppersmith_univariate
from sage.all import Integer, random_prime
# 5,7,11


def setup_rsa_and_cipher(known_prefix, unknown_length, e, m):
    N = random_prime(2**150) * random_prime(2**150)
    message = Integer(m, base=35)
    c = message**e % N
    return N, c, known_prefix, unknown_length, e


def test_coppersmith_univariate_e_3_1():
    N, c, known_prefix, unknown_length, e = setup_rsa_and_cipher(
        'thepasswordfortodayis000000000', 9, 3, 'thepasswordfortodayisswordfish')
    recovered_message = coppersmith_univariate.coppersmith_univariate(
        N, c, known_prefix, unknown_length, e)


def test_coppersmith_univariate_e_3_2():
    N, c, known_prefix, unknown_length, e = setup_rsa_and_cipher(
        '000000000aaa', 9, 3, 'bbbbbbbbbaaa')
    recovered_message = coppersmith_univariate.coppersmith_univariate(
        N, c, known_prefix, unknown_length, e)


def test_coppersmith_univariate_e_3_3():
    N, c, known_prefix, unknown_length, e = setup_rsa_and_cipher(
        'bb000000000aaa', 9, 3, 'bbcccccccccaaa')
    recovered_message = coppersmith_univariate.coppersmith_univariate(
        N, c, known_prefix, unknown_length, e)


# def test_coppersmith_univariate_e_3_2():
#     N, c, known_prefix, unknown_length, e = setup_rsa_and_cipher(
#         'thepasswordfortodayiss00000000', 8, 3)
#     recovered_message = coppersmith_univariate.coppersmith_univariate(
#         N, c, known_prefix, unknown_length, e)


# def test_coppersmith_univariate_e_3_3():
#     N, c, known_prefix, unknown_length, e = setup_rsa_and_cipher(
#         'thepasswordfortodayissw0000000', 7, 3)
#     recovered_message = coppersmith_univariate.coppersmith_univariate(
#         N, c, known_prefix, unknown_length, e)


# def test_coppersmith_univariate_e_3_4():
#     N, c, known_prefix, unknown_length, e = setup_rsa_and_cipher(
#         'thepasswordfortodayisswo000000', 6, 3)
#     recovered_message = coppersmith_univariate.coppersmith_univariate(
#         N, c, known_prefix, unknown_length, e)


# def test_coppersmith_univariate_e_3_5():
#     N, c, known_prefix, unknown_length, e = setup_rsa_and_cipher(
#         'thepasswordfortodayisswor00000', 5, 3)
#     recovered_message = coppersmith_univariate.coppersmith_univariate(
#         N, c, known_prefix, unknown_length, e)


# def test_coppersmith_univariate_e_3_6():
#     N, c, known_prefix, unknown_length, e = setup_rsa_and_cipher(
#         'thepasswordfortodayissword0000', 4, 3)
#     recovered_message = coppersmith_univariate.coppersmith_univariate(
#         N, c, known_prefix, unknown_length, e)


if __name__ == '__main__':
    pytest.main()
