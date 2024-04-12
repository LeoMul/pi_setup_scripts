import json 
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-j', '--json',  help='path of json')
args = parser.parse_args()



class Input:
    def __init__(self,name='defaultname',age='defaultage'):
        self.name = name
        self.age = age 

if not args.json:
    input_default = Input()
    default = json.dumps(input_default.__dict__,indent=1)
    print(default) 

else: 
    with open(args.json, 'r') as j:
        contents = json.loads(j.read())

    input = Input(**contents)

    print(input.age)
