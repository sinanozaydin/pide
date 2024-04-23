import pide
import sys

p_obj = pide.pide()
methods = [method for method in dir(p_obj) if callable(getattr(p_obj, method)) and not method.startswith(("__","_"))]
methods_set = [method for method in dir(p_obj) if callable(getattr(p_obj, method)) and method.startswith("set")]

print(methods_set[1])
print(p_obj.set_bulk_water.__doc__)
# Get the docstrings of methods



sys.exit()
print(method_docs)
# Print the docstrings
for name, docstring in method_docs.items():
	print(f"Method: {name}, Docstring: {docstring}")

# print(methods_set)