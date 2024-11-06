#!/usr/bin/env python
# libraries
## batteries
import os
import sys
import logging
import argparse
## 3rd party
from dotenv import load_dotenv
## package
from recharge_check.cli.employees import employees_main, employees_parser
from recharge_check.cli.inventory import inventory_main, inventory_parser
from recharge_check.cli.vivarium import vivarium_main, vivarium_parser
from recharge_check.utils import CustomFormatter

# logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)

# functions
def arg_parse(args=None) -> dict:
    """
    Parse command line arguments.
    """
    desc = "Recharge Check"
    epi = """DESCRIPTION:
    A tool for auto-checking issues with the Recharge system.
    """
    # check for OP
    if os.getenv("OPENAI_API_KEY") is None:
        raise ValueError("OPENAI_API_KEY not found in environment")
    
    # main parser
    parser = argparse.ArgumentParser(
        description=desc,
        epilog=epi,
        formatter_class=CustomFormatter
    )
    subparsers = parser.add_subparsers(dest="command", help="Subcommands")
    
    # subparsers
    ## employees
    employees_parser(subparsers)
    ## inventory
    inventory_parser(subparsers)
    ## vivarium
    vivarium_parser(subparsers)
    
    # parsing args
    return parser.parse_args()

def main():
    # load environment variables
    load_dotenv()
    # parsing args
    args = arg_parse()
    
    # which subcommand
    if args.command == "dataset-search":
        dataset_search_main(args)
    elif args.command == "inventory":
        inventory_main(args)
    elif args.command == "vivarium":
        vivarium_main(args)
    else:
        print("No command specified. Exiting ...")
        sys.exit(0)
    
if __name__ == "__main__":
    main()