'''
run_subprocess() function that unifies how
external tools are called within meningotype
'''

import subprocess
import sys

def run_subprocess(command):
    '''
    A function that takes a command line argument as a list
    and runs it, capturing any errors that the called command
    may produce.
    '''

    try:
        result = subprocess.run(command, capture_output=True, text=True, check=True)

        # print(f"Process '{' '.join(command)}' finished correctly.")

        cmd_output = result.stdout.strip()
        # print(f'Process output:\n{cmd_output}')
        return cmd_output

    except subprocess.CalledProcessError as err:
        cmd = f"'{' '.join(err.cmd)}'"
        print(f'ERROR running command {cmd}.\nCommand failed with exit code {err.returncode}:')
        print(f'{err.stderr.strip()}')
        sys.exit(err.returncode)

def main():
    '''
    Read commands from command line when called as stand-alone script
    '''

    command = sys.argv[1:]
    run_subprocess(command)

if __name__ == "__main__":
    main()
