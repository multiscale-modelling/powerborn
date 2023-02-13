import re 
import sys

def instrumentLine(line, accesses, line_number, lower_bound, upper_bound, array_name):
    # test if there is a declaration of a variable
    # type name; this one is not interesting since it has no accesses
    # type name = ...;
    # type name();
    if line.lstrip().startswith("__local"):
        return line + '\n'
    if line.find("//") < line.find(array_name) and line.find("//") >= 0:
        # we have a comment
        print line 
        return line + "\n"
    has_declaration = "=" in line and len(line.split("=")[0].split()) == 2
    is_if_statement = line.lstrip().startswith("if(") or line.lstrip().startswith("if (")
    if has_declaration:
        has_declaration = not array_name in line.split("=")[0].split()[0] 
    new_line = ""
    if has_declaration:
        new_line += line.split("=")[0] + ";\n" # the declaration
        new_line += "{\n" + " ".join(line.split()[1:]) + "\n" # the initialization
    elif not is_if_statement:
        new_line += "{\n" + line + "\n"
    elif is_if_statement:
        new_line += "{\n"
    for access, addr in accesses:
        if addr:
            new_line += '''if((%s < %s) || (%s >= %s)) printf((__constant char *) "Seg fault in line %d: %s at %s", %s, %s, %s, get_group_id(0));''' % (
                addr, lower_bound, addr, upper_bound, line_number, access, "index %d, lower %d, upper %d, group_id %d\\n", addr, lower_bound, upper_bound) + "\n"
    new_line += "}\n"
    if is_if_statement:
        new_line += line + "\n"
    return new_line

def instrumentKernel(array_name, src, lower_bound="0", upper_bound="0"):
    src_lines = src.splitlines()
    new_src = ""
    is_comment = False
    for i, line in enumerate(src_lines):
        if not is_comment and line.lstrip().startswith("/*"):
            is_comment = True           
        
        accesses = re.findall("(%s\[(.*?)\])" % array_name, line, re.DOTALL)
        
        if accesses and not is_comment:
            new_src += instrumentLine(line, accesses, i, lower_bound, upper_bound, array_name)
        else:
            if is_comment:
                print line
            new_src += line + "\n"
        
        if is_comment and line.rstrip().endswith("*/"):
            is_comment = False
    
    return new_src
    
    
if __name__ == "__main__":
    src_file = sys.argv[1]
    new_src_file = sys.argv[2]
    array_name = sys.argv[3]
    lower_bound = sys.argv[4]
    upper_bound = sys.argv[5]
    
    with open(src_file) as fin:
        src = fin.read()

    new_src = instrumentKernel(array_name, src, lower_bound, upper_bound)
    with open(new_src_file, "w") as fout:
        fout.write(new_src)
