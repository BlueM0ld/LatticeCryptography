from utility import delete_folder,update_filename
update_filename("../src") #pre execution
############################################
from graph_plotting import combine_gso_norms, generate_graphs
from coppersmith_univariate import coppersmith_univariate, setup_rsa_and_cipher
from coron_direct_bivariate import coron_attack, setup
from hastad import hastads_attack_lattice, set_up_hastad
import os

valid_reductions = {"LLL", "BKZ", "BKZ40", "BKZ60"}
attack_options = [
        "coppersmith_univariate",
        "coron_direct_bivariate",
        "hastad",
        "weiner"
]


persist_state = False # this is a flag for persisting files

def check_choice(choice):
    if choice not in ['i', 'r']:
        print("Invalid choice! Please enter 'r' for random or 'i' for input.")
        sys.exit(1)  # Exit gracefully!

def check_reduction(reduction):
    if reduction not in valid_reductions:
        print(f"Invalid reduction type '{reduction}'. Must be one of {valid_reductions}")
        sys.exit(1)  # Exit gracefully!

def check_graph_choice(choice):
    if choice>len(attack_options) and choice <0:
        print(f"Invalid choice! Muse be 1 - {len(attack_options)}")
        sys.exit(1)  # Exit gracefully!


def persist_results():
    global persist_state 
    # Toggle the state
    persist_state = not persist_state
    
    print("-"*50)
    print("This will persist results between runs, if this is not what you wanted, select this option again")

    print("Currently the persist option is:", persist_state)


def run_generate_graph():


    print("Which attack you want to generate graphs for:")
    print("1. Coppersmith's Univariate Attack")
    print("2. Coron's Direct Bivariate Attack")
    print("3. Hastad's Attack")
    print("4. Weiner's Attack")

    try:
        choice = int(input("Enter the number... "))
        check_graph_choice(choice)

        reduction = input("What type of reduction (LLL, BKZ, BKZ40, BKZ60): ").strip().upper()
        check_reduction(reduction)

        choice_g = int(input("Generate individual graphs or combine all? (Enter '1' for individuals or '2' for combined): "))
        filename = input("Enter file name (default will look for filename - reduced_gso_norms): ").strip()

        if choice_g == 1:
            generate_graphs(attack_options[choice-1], reduction, filename)
        elif choice_g == 2:
            combine_gso_norms(attack_options[choice-1], reduction, filename)
        else:
            print("Invalid choice. Please select a valid attack.")

    except ValueError:
        print("Invalid choice! Exiting....")
        return



def run_hastads_attack():
    #Not necessary to have random input choice as they are the same
    try:
        message = input("Enter the message: ")
        n = int(input("Enter the number of times to run the attack: "))
        bits = int(input("Enter number of bits: "))

        reduction = input("What type of reduction (LLL, BKZ, BKZ40, BKZ60): ").strip().upper()
        check_reduction(reduction)

        e = int(input("Enter e: "))
        k = int(input("Enter k (Number of individuals): "))

        for i in range(n):
            print(f"\n Run {i + 1}/{n}")

            ciphertexts, N = set_up_hastad(message, bits, e, k)
            try:
                recovered_message = hastads_attack_lattice(ciphertexts, N, e, reduction)
                print(f"recovered message with {k} recipients: {recovered_message}")
            except ValueError as ve:
                print(f"Failed with {k} recipients: {ve}")

    except ValueError:
        print("Invalid choice! Exiting....")
        return


def run_coron_bivariate():

    #Not necessary to have random input choice as they are the same
    try:
        n = int(input("Enter the number of times to run the attack: "))
        bits = int(input("Enter number of bits: "))

        reduction = input("What type of reduction (LLL, BKZ, BKZ40, BKZ60): ").strip().upper()
        check_reduction(reduction)

        k = int(input("Enter k: "))
        for i in range(n):
            print(f"\n Run {i + 1}/{n}")
            p, q, N, mask = setup(bits)
            coron_attack(mask, N,p,q , reduction, k=k)

    except ValueError:
        print("Invalid choice! Exiting....")
        return


def run_coppersmith_univariate():


    choice = input("Generate a random instance or input your own? (Enter 'r' for random or 'i' for input): ").strip().lower()
    check_choice(choice)

    try:
        n = int(input("Enter the number of times to run the attack: "))
        
        reduction = input("What type of reduction (LLL, BKZ, BKZ40, BKZ60): ").strip().upper()
        check_reduction(reduction)
        e = int(input("Enter e: "))
        bit_length = int(input("Enter number of bits: "))


        if choice == 'i':
            message = input("Enter the message: ")
            known_length = int(input("Enter the length of the known prefix: "))
            
            assert known_length <= len(message), "Error - length of the known prefix cannot exceed the length of the message!"


            known_prefix = message[:known_length]
            unknown_length = len(message) - known_length

            for i in range(n):
                print(f"\n Run {i + 1}/{n}")

                N, c, known_prefix, unknown_length, e = setup_rsa_and_cipher(known_prefix, unknown_length, e, message, bit_length)

                try:
                    recovered_message = coppersmith_univariate(N, c, known_prefix, unknown_length, e, reduction)
                    print(f"Recovered message with known prefix: {known_prefix}")
                    print(f"Was the recovered message correct? : {recovered_message == message}")
                except Exception as ex:
                    print(f"Failed with known prefix length {len(known_prefix)}: {ex}")

        elif choice == 'r':
            message_length = int(input("Enter the total length of the random message: "))
            known_length = int(input("Enter the length of the known prefix: "))

            assert known_length <= message_length, "Error - length of the known prefix cannot exceed the length of the message!"


            for i in range(n):
                print(f"\n Run {i + 1}/{n}")

                message = os.urandom(message_length).hex()
                known_prefix = message[:known_length]
                unknown_length = len(message) - known_length

                N, c, known_prefix, unknown_length, e = setup_rsa_and_cipher(known_prefix, unknown_length, e, message, bit_length)

                try:
                    recovered_message = coppersmith_univariate(N, c, known_prefix, unknown_length, e, reduction)
                    print(f"Recovered message with known prefix: {known_prefix}")
                    print(f"Was the recovered message correct? : {recovered_message == message}")
                except Exception as ex:
                    print(f"Failed with known prefix length {len(known_prefix)}: {ex}")
        else:
            print("Invalid choice! Exiting....")
            return
    except ValueError:
        print("Invalid choice! Exiting....")
        return
       

def main():
    """
    Main interface for the operations that can be executed from the script
    """

    while True:
        print("-"*50)
        print("Options available:")
        print("1. Coppersmith Univariate Attack")
        print("2. Coron's Bivariate Attack")
        print("3. Hastad's Attack")
        print("4. Weiner's Attack")
        print("5. Generate Graphs")
        print("6. Delete Results")
        print("7. Persist Results")
        print("8. Exit")



        choice = input("Enter a number to select and run: ")

        if persist_state is True:
            update_filename("../result")

        if choice == '1':
            run_coppersmith_univariate()
        elif choice == '2':
            run_coron_bivariate()
        elif choice == '3':
            run_hastads_attack()
        elif choice == '4':
            run_weiners_attack()
        elif choice == '5':
            run_generate_graph()
        elif choice == '6':
            delete_folder()
        elif choice == '7':
            persist_results()
        elif choice == '8':
            sys.exit("Exiting...bye!")
        else:
            print("Invalid choice. Please select a valid attack.")

if __name__ == "__main__":

    main()
