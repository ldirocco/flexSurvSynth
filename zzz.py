import textwrap

def package_startup_message():
    info = "Find out more at https://www.synthpop.org.uk/"
    wrapped_info = textwrap.fill(info, width=80)
    print(wrapped_info)

# Call the function when the module is imported
package_startup_message()
