from dolfin import *

for x in  parameters.keys():
    print x

print '\nform_compiler parameters'

for x in parameters['form_compiler'].keys():
    print x