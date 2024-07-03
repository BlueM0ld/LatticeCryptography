import sys
import os

# Add the path to the src directory to sys.path
src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(src_path)

from sage.all import Integer, random_prime
from coppersmith_univariate import coppersmith_univariate


def setup_rsa_and_cipher(known_prefix, unknown_length, e, m):
    N = random_prime(2**150) * random_prime(2**150)
    message = Integer(m, base=35) #  0-9 - a-z
    c = message**e % N
    return N, c, known_prefix, unknown_length, e


def experiment_coppersmith__e_3_1():
    N, c, known_prefix, unknown_length, e = setup_rsa_and_cipher(
        'thepasswordfortodayis000000000', 9, 3, 'thepasswordfortodayisswordfish')
    recovered_message = coppersmith_univariate(
        N, c, known_prefix, unknown_length, e)


def experiment_coppersmith__e_3_2():
     N, c, known_prefix, unknown_length, e = setup_rsa_and_cipher(
         '000000000aaa', 9, 3, 'bbbbbbbbbaaa')
     recovered_message =coppersmith_univariate(
         N, c, known_prefix, unknown_length, e)


def experiment_coppersmith__e_3_3():
     N, c, known_prefix, unknown_length, e = setup_rsa_and_cipher(
         'bb000000000aaa', 9, 3, 'bbcccccccccaaa')
     recovered_message =coppersmith_univariate(
         N, c, known_prefix, unknown_length, e)


def experiment_coppersmith__e_3_4():
     N, c, known_prefix, unknown_length, e = setup_rsa_and_cipher(
         'hello000000000', 10, 3, 'helloearthlings1')
     recovered_message =coppersmith_univariate(
         N, c, known_prefix, unknown_length, e)


def experiment_coppersmith__e_3_5():
     N, c, known_prefix, unknown_length, e = setup_rsa_and_cipher(
         'c000000didi0', 7, 3, 'charliedidit')
     recovered_message =coppersmith_univariate(
         N, c, known_prefix, unknown_length, e)


def experiment_coppersmith__e_3_6():
     N, c, known_prefix, unknown_length, e = setup_rsa_and_cipher(
         'ali0ea0db0b', 3, 3, 'aliceandbob')
     recovered_message =coppersmith_univariate(
         N, c, known_prefix, unknown_length, e)


def experiment_coppersmith__e_3_7():
     N, c, known_prefix, unknown_length, e = setup_rsa_and_cipher(
         '123056789000054321', 5, 3, '123456789987654321')
     recovered_message =coppersmith_univariate(
         N, c, known_prefix, unknown_length, e)


def experiment_coppersmith__e_3_8():
     N, c, known_prefix, unknown_length, e = setup_rsa_and_cipher(
         'ihate00000e', 5, 3, 'ihatecheese')
     recovered_message =coppersmith_univariate(
         N, c, known_prefix, unknown_length, e)


def experiment_coppersmith__e_5_1():
     N, c, known_prefix, unknown_length, e = setup_rsa_and_cipher(
         'thepasswordfortodayis000000000', 9, 5, 'thepasswordfortodayisswordfish')
     recovered_message =coppersmith_univariate(
         N, c, known_prefix, unknown_length, e)


def experiment_coppersmith__e_5_2():
     N, c, known_prefix, unknown_length, e = setup_rsa_and_cipher(
         'thepasswordfortodayisswo000000', 6, 5, 'thepasswordfortodayisswordfish')
     recovered_message =coppersmith_univariate(
         N, c, known_prefix, unknown_length, e)


def experiment_coppersmith__e_5_3():
     N, c, known_prefix, unknown_length, e = setup_rsa_and_cipher(
         'thepasswordfortodayisswordf000', 3, 5, 'thepasswordfortodayisswordfish')
     recovered_message =coppersmith_univariate(
         N, c, known_prefix, unknown_length, e)


def experiment_coppersmith__e_5_4():
     N, c, known_prefix, unknown_length, e = setup_rsa_and_cipher(
         'b00c0cc0c0ca0', 6, 5, 'bbcccccccccaaa')
     recovered_message =coppersmith_univariate(
         N, c, known_prefix, unknown_length, e)

def experiment_coppersmith__e_5_5():
     N, c, known_prefix, unknown_length, e = setup_rsa_and_cipher(
         'bbccccccccc000', 3, 5, 'bbcccccccccaaa')
     recovered_message =coppersmith_univariate(
         N, c, known_prefix, unknown_length, e)

if __name__ == '__main__':
    experiment_coppersmith__e_3_1()
    experiment_coppersmith__e_3_2()
    experiment_coppersmith__e_3_3()
    experiment_coppersmith__e_3_4()
    experiment_coppersmith__e_3_5()
    experiment_coppersmith__e_3_6()
    experiment_coppersmith__e_3_7()
    experiment_coppersmith__e_3_8()
    experiment_coppersmith__e_5_1()
    experiment_coppersmith__e_5_2()
    experiment_coppersmith__e_5_3()
    experiment_coppersmith__e_5_4()
    experiment_coppersmith__e_5_5()
#    experiment_coppersmith__e_5_6()
#    experiment_coppersmith__e_5_7()
#    experiment_coppersmith__e_5_8()
#    experiment_coppersmith__e_5_9()
#    experiment_coppersmith__e_5_10()
#    experiment_coppersmith__e_5_11()
#    experiment_coppersmith__e_5_12()
#    experiment_coppersmith__e_5_13()
