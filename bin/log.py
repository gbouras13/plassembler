import logging
import sys

def write_to_log(s, logger):
           while True:
                output = s.readline().decode()
                if output:
                    logger.log(logging.INFO, output)
                else:
                    break


def write_message(message, logger):
    print(message)
    logger.info(message)

def print_and_exit(message, logger):
    logger.info(message)
    sys.exit(message)