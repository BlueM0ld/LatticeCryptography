from sage.all import Integer, log, matrix, random_prime, RR, ZZ
from fpylll import GSO, IntegerMatrix

def coppersmith_univariate():
    print("Executing Coppersmith's Univariate Attack...")
    # TODO: Implement coppersmith_univariate attack
    pass

def coppersmith_bivariate():
    print("Executing Coppersmith's Bivariate Attack...")
    # TODO: Implement coppersmith_bivariate attack
    pass

def coppersmith_multivariate():
    print("Executing Coppersmith's Multivariate Attack...")
    # TODO: Implement coppersmith_multivariate attack
    pass

if __name__ == '__main__':
    choice = input("Would you like to generate random primes or input your own? (random/input): ").strip().lower()
    
    if choice == 'input':
        p = Integer(input('Enter prime number p: ').strip())
        q = Integer(input('Enter prime number q: ').strip())
    else:
        print('Generating primes')
        p = random_prime(2^1024)
        q = random_prime(2^1024)

    print('Calculating modulus N')
    N = p * q
    print(f'N = {N}')
    
    method_choice = input("Which Coppersmith method would you like to call? (univariate/bivariate/multivariate): ").strip().lower()
    
    if method_choice == 'univariate':
        coppersmith_univariate()
    elif method_choice == 'bivariate':
        coppersmith_bivariate()
    elif method_choice == 'multivariate':
        coppersmith_multivariate()
    else:
        print("Invalid choice. Please select from 'univariate', 'bivariate', or 'multivariate'.")
