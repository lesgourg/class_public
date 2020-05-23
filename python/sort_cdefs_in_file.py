file = 'classy.pyx'
with open(file, 'r') as fid:
    lines = fid.readlines()
    cdef_block_found = False
    current_block_lines = []
    for line_index, line in enumerate(lines):
        if not cdef_block_found:
            cdef_position = line.find('cdef:')
            if cdef_position:
                cdef_block_found = True
                indentation = cdef_position + 4
        else:
            # We are in a block
            if len(line) <= indentation or len(line[:indentation].strip()) != 0 or line[indentation] == ' ':
                # The block has ended
                cdef_block_found = False
                # Sort the block
                block_length = len(current_block_lines)
                lines[line_index - block_length:line_index] = sorted(current_block_lines)
                current_block_lines = []
            else:
                current_block_lines.append(line)

with open(file, 'w') as fid:
    fid.writelines(lines)