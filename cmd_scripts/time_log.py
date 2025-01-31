import os
import pandas as pd
from datetime import datetime
import argparse

parser = argparse.ArgumentParser(
    description='Add a new line to the time log file...')
parser.add_argument('--message', help='message to add to the time log')
parser.add_argument('--note_path', help='path to the time log file')

args = parser.parse_args()

if args.message is None or args.note_path is None:
    print("Both --message and --note_path  are required.")
    exit(1)

print("\n\nmessage: ", args.message, "\npath: ",
      args.note_path)
log_file = args.note_path
if os.path.exists(log_file):
    df = pd.read_csv(log_file)
else:
    df = pd.DataFrame(columns=["datetime", "comment"])

row = {"datetime": datetime.today(), "comment": args.message}
df = pd.concat([df, pd.DataFrame(row, index=[0])], ignore_index=True)
df.to_csv(log_file, index=False)
