from m2r2 import convert


def build_quick_start():

    replace = {'# pyremap': '# Quick Start'}

    skip = [('## Documentation', '## Installation')]
    outContent = ''
    skipMode = False
    with open('../README.md', 'r') as inFile:
        for line in inFile.readlines():
            for skipStart, skipEnd in skip:
                if not skipMode and skipStart in line:
                    skipMode = True
                if skipMode and skipEnd in line:
                    skipMode = False
            if not skipMode:
                for replaceString in replace:
                    if replaceString in line:
                        line = replace[replaceString]
                        break
                outContent = outContent + line

    outContent = convert(outContent)

    with open('quick_start.rst', 'w') as outFile:
        outFile.write('.. _quick_start:\n\n')
        outFile.write(outContent)
