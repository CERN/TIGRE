import os
def run():
    dirname = os.path.dirname(__file__)
    dirname = os.path.join(dirname,'launch.sh')
    os.system('bash '+dirname)
