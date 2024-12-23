# Helper function to go from str -> int or float

def parser(file):
    def to_numeric(in_str):
        try:
            return int(in_str)
        except:
            return float(in_str)

    file = file
    parsed = {}
    with open(file, 'r') as f:
        is_parsing = False
        is_2d = False # 2D table
        is_3d = False # "3D" table, i.e. a table w/ 2 indices
        is_simple = False  # is_simple -> := on first line
        variable_name = ''
        variable_value = None
        parsing_value = None

        for line in f.readlines():
            # Split the line into tokens
            tokens = list(filter(None, line.strip('\n').replace('\t', ' ').split(' ')))

            # Skip over empty lines
            if len(tokens) == 0:
                continue

            # Remove comments (anything after and including the #
            for token_idx in range(len(tokens)):
                if '#' in tokens[token_idx]:
                    tokens = tokens[:token_idx]
                    break

            if tokens[0] == 'param':
                # New parameter!
                variable_name = tokens[1].split(':')[0]
                is_parsing = True
                # Is this a single-line assignment (i.e. do we need to parse the value too)?
                for token in tokens:
                    if ':=' in token:
                        # Value coming up (variables like TF:= and then continue on next line) -> need to set parsing_value so we actually process those lines that follow, without first hitting a header
                        parsing_value = []
                    elif ';' in token:
                        value = token.strip(';')
                        parsed[variable_name] = to_numeric(value)
                        is_parsing = False
            elif is_parsing and ':=' in tokens:  # This is a header line
                # Multiple indices -> push to first index; only support 2 indexing variables at most for now
                if parsing_value is not None and parsing_value != []:
                    if variable_value is None:
                        variable_value = []
                    # Add a copy!
                    variable_value.append(parsing_value.copy())

                # Always reset for fresh data
                parsing_value = []
            elif [';'] == tokens:
                # End of parse
                # Add parsing_value to variable_value
                if parsing_value is not None:
                    if variable_value is None:
                        variable_value = parsing_value
                    else:
                        variable_value.append(parsing_value)
                # Add it to the parsed dict
                parsed[variable_name] = variable_value

                # Reset variables
                parsing_value = None
                variable_value = None
                variable_name = None

            elif parsing_value is not None:
                # Strip first column (non-numerical generally), and add the rest
                # TODO: turn into int or float instead?
                line_values = list(to_numeric(token) for token in tokens[1:])
                if len(line_values) == 1:
                    parsing_value.append(line_values[0]) # "Unwrap" the line's values so "1D" tables don't have [[1], [2], ...] but can be indexed with 1 variable (i.e. [1, 2, ...]) instead
                else:
                    parsing_value.append(line_values)
            else:
                pass
                #print(f"Unable to process token string: {tokens}")

            # TODO: more parsing
            #print(tokens)

    #print(parsed.keys())
    #print(parsed)
    # print('===')
    # print(parsed['E'])
    return parsed