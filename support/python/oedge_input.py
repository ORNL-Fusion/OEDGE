try:
    # framework is running
    from .startup_choice import *
except ImportError as _excp:
    # class is imported by itself
    if (
        'attempted relative import with no known parent package' in str(_excp)
        or 'No module named \'omfit_classes\'' in str(_excp)
        or "No module named '__main__.startup_choice'" in str(_excp)
    ):
        from startup_choice import *
    else:
        raise

from omfit_classes.omfit_ascii import OMFITascii
from omfit_classes.omfit_nc import OMFITnc
from omfit_classes import namelist

from collections import OrderedDict

__all__ = ['OMFIToedgeInput', 'OMFIToedgeRun', 'OMFIToedgeNC']


def interpret(x):

    x = x.split('\'')
    for k in range(1, len(x), 1)[::2]:
        x[k] = '\'' + x[k] + '\''
    x = [_f for _f in x if _f]

    xx = []
    for k in x:
        if k.startswith('\''):
            xx.append(k)
        else:
            xx.extend(list(map(namelist.interpreter, [_f for _f in re.split(' |\t', k) if _f])))

    data = []
    inline_comment = ''

    for k in range(len(xx)):
        if not len(inline_comment) and is_numeric(xx[k]):
            data.append(xx[k])
        elif not len(inline_comment) and xx[k][0] == '\'':
            data.append(xx[k][1:-1])
        else:
            inline_comment += ' ' + str(xx[k])
    inline_comment = inline_comment.strip(' \'')

    return data, inline_comment


class OMFIToedgeInput(SortedDict, OMFITascii):
    """OEDGE input data file"""

    def __init__(self, filename, codefilename=None, casename=None, caseloc=None, root=None, **kw):
        OMFITascii.__init__(self, filename, **kw)
        SortedDict.__init__(self)

        self.dynaLoad = True

        printd('Init OEDGE input file: ', casename, caseloc)

        self.casename = casename
        self.codefilename = codefilename
        self.caseloc = caseloc

        self.tagseries = []
        self.subs = []
        
        # set casename and treelocation of the object
        self.setloc(casename, caseloc, root)

    def __call__(self, *args):
        # Allow the case object to return the string representation of its tree elements
        # by checking case is eval(case()) in the calling routine - can verify object is loaded where expected
        if len(args) == 0:
            return self.caseloc
        elif self.caseloc is not None:
            tmp = self.caseloc
            for arg in args:
                tmp = tmp + "['" + arg + "']"
            return tmp

    def setloc(self, casename=None, caseloc=None, root=None):
        # Rotuine will assign the casename and location in the tree either to passed in values
        # or it will attempt to determine them based on the filename and treelocation return values
        # print('filename:',self.filename)
        # print('location:',treeLocation(self))
        # print('inputs:',casename,caseloc)
        # print('current:',self.casename,self.caseloc)

        # Allow overwrite of case name  ... do I need to update filename and write new file (?)
        # leave as a separate case in the IF statement so this can be addressed later
        if self.casename is None and casename is None:
            self.casename, filetype = os.path.splitext(self.filename.split(os.sep)[-1])
        elif self.casename is None and casename is not None:
            self.casename = casename
        elif casename is not None:
            self.casename = casename

        # Case is not in the tree when first created on module load .. object is added to the
        # tree after it is created (of course) so treeLocation is not usable until later.

        if self.caseloc is None and caseloc is None:
            #  check correct location in tree first ... then use treeLocation if OMFIT not available and root not specified
            if self.casename is not None:
                try:
                    # try is needed because dictionary entry may not be accessible at module load
                    # If root is not specified then try using the absolute path
                    if root is None:
                        if self.casename in OMFIT['OEDGE']['INPUTS']['CASES']:
                            # check to see if the object is loaded into the tree in the expected location
                            if self is OMFIT['OEDGE']['INPUTS']['CASES'][self.casename]:
                                self.caseloc = "OMFIT['OEDGE']['INPUTS']['CASES'][%s]" % repr(self.casename)
                            else:
                                self.caseloc = None
                    else:
                        if self.casename in root['INPUTS']['CASES']:
                            # check to see if the object is loaded into the tree in the expected location
                            if self is root['INPUTS']['CASES'][self.casename]:
                                self.caseloc = "root['INPUTS']['CASES'][%s]" % repr(self.casename)
                            else:
                                self.caseloc = None

                except KeyError:
                    self.caseloc = treeLocation(self)[-1]
                    # when the object is not loaded in the tree then treeLocation returns the null string
                    # replace with None
                    if self.caseloc == '':
                        self.caseloc = None

        elif self.caseloc is None and caseloc is not None:
            self.caseloc = caseloc
        elif caseloc is not None:
            self.caseloc = caseloc

    @dynaLoad
    def load(self):

        # empty the dictionary in case this is a reload
        self.clear()

        with open(self.filename, 'r') as f:
            lines = f.read()

        # handle inline comments with $
        lines = re.sub(r'(?<!\n)\$', '\n$', lines)
        lines = lines.split('\n')

        # split sections
        c = 0
        inhibit = False
        item_cnt = 0
        for line0 in lines:
            line = line0.rstrip()

            # stop parsing when find a line starting with c
            if line.startswith('c'):
                inhibit = True
                self['__comment_%06d__' % c] = line0
                c += 1
            # new entry with Z'{ style
            elif not inhibit and line.startswith(r'\'{'):
                item = re.sub(r'.*(\{.*\}).*', r'\1', line.split("'")[1])
                self[item] = SortedDict()
                self[item]['__raw__'] = [line]
            # new entry with '* or '+ style
            elif not inhibit and re.match(r'^\'[\*\+][a-zA-Z0-9][0-9][0-9]', line):
                item = line.split("'")[1][:4]
                item_cnt += 1
                # Item C14 is repeated twice in some input files - the second should be C15
                if item in self:
                    # found a repeated item
                    # print("Repeat item found:",item)
                    # if the repeat item is C14 replace with C15
                    if item == '+C14':
                        item == '+C15'
                self[item] = SortedDict()
                self[item]['__raw__'] = [line]
            # continuation of entry
            elif not inhibit and re.match(r'^ *[\'0-9\.\-].*', line):
                # Think this is where the other mislabeled lines go
                #
                # Kludgy fix to a bug in older input files
                #
                #if 'Ring Range :-' in line:
                #    # insert tag at beginning of line
                #    item = '+F07'
                #    line = "'" + item + line[1:]
                #    self[item] = SortedDict()
                #    self[item]['__raw__'] = [line]
                #elif 'J0 & J1    :-' in line:
                #    item = '+F08'
                #    line = "'" + item + line[1:]
                #    self[item] = SortedDict()
                #    self[item]['__raw__'] = [line]
                #elif 'VMF0,1,2   :-' in line:
                #    item = '+F09'
                #    line = "'" + item + line[1:]
                #    self[item] = SortedDict()
                #    self[item]['__raw__'] = [line]
                #else:
                self[item]['__raw__'].append(line)
            # comment
            else:
                self['__comment_%06d__' % c] = line0
                c += 1

        print("Total tags in input file:",item_cnt)
                
        # interpret sections
        for item in list(self.keys()):
            if item.startswith('__comment_'):
                continue

            # initialize the tag field for all non-comment input items
            self[item]['tag'] = ''

            raw = self[item]['__raw__']
            raw0 = raw[0].strip().split('\'')
            # special {...}
            if '{' in item:
                self[item]['comment'] = raw0[1].split('}')[1].strip()
                self[item]['data'] = []
                self[item]['tag'] = item
                self[item]['inline_comment'] = ''

                try:
                    # scalar
                    if not raw0[2]:
                        raise
                    self[item]['data'], inline_comment = interpret(raw0[2])
                    # scalars as scalars
                    if len(self[item]['data']) == 1:
                        self[item]['data'] = self[item]['data'][0]

                except Exception:
                    # array
                    for line in raw[1:]:
                        data, inline_comment = interpret(line)
                        self[item]['data'].append(data)
            # standard
            else:
                self[item]['comment'] = raw0[1][5:].strip()
                self[item]['tag'] = item[1:]

                # scalar(ish)
                if len(raw) == 1:
                    tmp = interpret(raw[0])
                    self[item]['data'] = tmp[0][1:]
                    self[item]['inline_comment'] = tmp[-1]
                    # scalars as scalars
                    if len(self[item]['data']) == 1:
                        self[item]['data'] = self[item]['data'][0]

                # array
                else:
                    self[item]['inline_comment'] = '\''.join(raw0[2:])
                    self[item]['inline_comment'] = self[item]['inline_comment'].strip(' \'')
                    for ndummies, line in enumerate(raw[1:]):
                        if "' '" not in line:
                            break
                    header = re.sub(r'(\'.*\'\s*)([0-9]+)(\s.*)', r'\1', raw[1 + ndummies] + ' ').strip('\' ')
                    n = int(re.sub(r'(\'.*\'\s*)([0-9]+)(\s.*)', r'\2', raw[1 + ndummies] + ' ').strip('\' '))
                    header_comment = re.sub(r'(\'.*\'\s*)([0-9]+)(\s.*)', r'\3', raw[1 + ndummies] + ' ').strip()
                    data = []
                    for k in range(len(raw[2 + ndummies :])):
                        data.append(interpret(raw[2 + ndummies + k])[0])
                    dummies = raw[1 : ndummies + 1]
                    self[item]['dummies'] = dummies
                    self[item]['header'] = header
                    self[item]['n'] = n
                    self[item]['header_comment'] = header_comment
                    self[item]['data'] = data



        #if (self.codefilename is not None):
        #    load_code(self,self.codefilename)
                
            

    def load_code(self,incode):

        # this loads supplemental data from the code file


        # if a supplemental file with code for the tags has been specified then read it in and process it.
        with open(incode, 'r') as f:
            lines = f.read()
            

        # handle inline comments with $
        lines = re.sub(r'(?<!\n)\$', '\n$', lines)
        lines = lines.split('\n')

        # split sections
        c = 0
        inhibit = False
        item = None
        item_cnt = 0
        for line0 in lines:
            line = line0.strip()

            # include comments in code file
            if line.startswith('!'):
                if item is not None:
                    if '__coderaw__' in self[item]:
                        self[item]['__coderaw__'].append(line) 

            #continue
            # new entry with Z'{ style
            elif not inhibit and re.match(r'^\'[\*\+][a-zA-Z0-9][0-9][0-9]', line):
                item = line.split("'")[1][:4]
                item_cnt += 1
                if item in self:
                    #print("Start:",item,':',line.strip())
                    self[item]['__coderaw__'] = [line.strip()]
                else:
                    print("READING CODE: TAG NOT IN INPUT: ",item)
            # continuation of entry
            elif not inhibit and re.match(r'^ *[\'0-9\.\-].*', line):
                if item in self:
                    #print("Add A:",item,':',line.strip())                    
                    self[item]['__coderaw__'].append(line.strip())
            else:
                if item in self:
                    #print("Add B:",item,':',line.strip())                    
                    self[item]['__coderaw__'].append(line.strip())

        print("Total tags with code:",item_cnt)
                    
        for item in list(self.keys()):        
            if item.startswith('__comment_'):
                continue
            # process the coderaw
            if '__coderaw__' in self[item]:
                #print('coderaw: ', item, ' coderaw: ', self[item]['__coderaw__'])
                self[item]['code_input'] = []
                self[item]['code_comments'] = []
                for code in list(self[item]['__coderaw__']):
                    if len(code) < 1:
                        continue
                    #print("Item:",item," code:",code[0:4].lower(),': ',code)
                    if code[0:4].lower() == 'call':
                        # arguments without brackets
                        args =  code.split('(',1)[1][:-1]
                        sub =   (code.split('(',1)[0]).split(' ',1)[1].lower()
                        self[item]['code'] = code
                        self[item]['sub'] = sub.strip()
                        self[item]['args'] = args
                        #print('Extracted:',item,'Sub:',self[item]['sub'],'Args:',self[item]['args'])
                    elif code[0:1] == '!':
                        self[item]['code_comments'].append(code)                            
                    else:
                        # collect everything except the code line and blank lines
                        self[item]['code_input'].append(code)
                # process subrotuine call arguments
                if 'sub' in self[item]: 
                    sub = self[item]['sub']
                    args = self[item]['args']

                    values = args.split(',')

                    if sub == 'rdi' or sub == 'rdr' or sub == 'rdq':

                        self[item]['vars']  = [values[0]]
                        self[item]['var_doc'] = [values[5]]

                    elif sub == 'rdi2' or sub == 'rdr2' or sub == 'rdq2':

                        self[item]['vars']  = [values[0],values[1]]
                        self[item]['var_doc'] = [values[6]]

                    elif sub == 'rdi3' or sub == 'rdr3' or sub == 'rdq3':

                        self[item]['vars']  = [values[0],values[1],values[2]]
                        self[item]['var_doc'] = [values[7]]

                    elif sub == 'rdc':
                        self[item]['vars']  = [values[0]]
                        self[item]['var_doc'] = [values[1]]
                        
                    elif sub == 'rdrar':
                        
                        self[item]['vars']  = [values[1]]
                        self[item]['var_doc'] = [values[6]]
                        self[item]['var_arr']  = [values[0]]
                        
                    elif sub == 'rdrarn':

                        self[item]['vars']  = [values[1]]
                        self[item]['var_doc'] = [values[9]]
                        self[item]['var_arr']  = [values[0]]
                        
                    elif sub == 'rdvmf':

                        self[item]['var_doc'] = [values[0]]
                        
                    else:
                        print("Args: Unknown subroutine found:",sub)


                # process sample input to obtain default values for arguments to subroutine calls
                if 'code_input' in self[item]:
                    if 'sub' in self[item]:
                        sub = self[item]['sub']

                        inline = self[item]['code_input'][0]
                    
                        if sub == 'rdi' or sub == 'rdr' or sub == 'rdq':
                            
                            defvals = inline.rsplit('\'',1)[1].split()
                            self[item]['defvals']  = defvals

                        elif sub == 'rdi2' or sub == 'rdr2' or sub == 'rdq2':

                            defvals = inline.rsplit('\'',1)[1].split()                        
                            self[item]['defvals']  = defvals

                        elif sub == 'rdi3' or sub == 'rdr3' or sub == 'rdq3':

                            defvals = inline.rsplit('\'',1)[1].split()
                            self[item]['defvals']  = defvals

                        elif sub == 'rdc':

                            # extract second string on line
                            defvals = inline[inline[0:len(inline)-1].rfind('\'',1):]
                            self[item]['defvals']  = [defvals]
                        
                        elif sub == 'rdrar' or sub == 'rdrarn':
                            
                            inline2 = self[item]['code_input'][1]
                            defvals = inline2.rsplit('\'',1)[1].split()                        
                            self[item]['defvals']  = defvals
                            print("ITEM:",item,":sub:",sub,":defvals:",self[item]['defvals'])
                        
                        elif sub == 'rdvmf':

                            self[item]['defvals'] = ['special']
                        
                        

                    else:
                        print("Input: Unknown subroutine found:",sub)
                    
                    
                        
        self.define_tag_series()


        

    def define_tag_series(self):

        # go through all of the tags - create a list of tag series - also check to see if any tags do not both have an input and code value
        tag_cnt = 0
        for item in list(self.keys()):
            # ignore comments
            if item.startswith('__comment_'):
                continue
            tag_cnt += 1
            tag = self[item]['tag']
            series = tag[0:1]
            if series not in self.tagseries:
                self.tagseries.append(series)

        #print("Total tags:",tag_cnt)
        #print("Tag series:",self.tagseries)

        sub_cnt = 0
        for item in list(self.keys()):
            # ignore comments
            if item.startswith('__comment_'):
                continue
            sub_cnt += 1
            tag = self[item]['tag']
            sub = self[item]['sub']
            if sub not in self.subs:
                self.subs.append(sub)

        print("Total subs:",sub_cnt)
        print("Subs      :",self.subs)
                

    def create_init_code(self,name):
        # generate the code to initialize the unstructured variables to default values - decide to use either values in the sample code or reference file input values

        #print("stub")

        # generate the code to initialize the inputs
        init_code = []
        # loop through tagseries
        item_cnt = 0
        for series in self.tagseries:
            
            series_cnt = 0
            init_code.append('')
            init_code.append('!')
            init_code.append('!  Initialization for TAG series %s'%(series))
            init_code.append('!')

            for item in self:
                # ignore comments
                if item.startswith('__comment_'):
                    continue

                # process tags matching the tag series - compose output file
                # contains sample input - commented original read code - converted read code

                # sample input -
                #
                tag = self[item]['tag']
                s = tag[0:1]
                
                if s == series:

                    # start with comments - tag heading - input example - any other documentation - then add init code

                    item_cnt += 1
                    # print out tag and all comments/documentation etc

                    sub = self[item]['sub']
                    
                    init_code.append('!')
                    init_code.append('! TAG:  %s    Initialization'%(item))
                    init_code.append('!')

                    # Add comments if any
                    if 'code_comments' in self[item]:
                        if len(self[item]['code_comments']) > 0:
                            init_code.append('! Comments from input')
                            init_code.append('!')
                            for comment in self[item]['code_comments']:
                               init_code.append(comment)
                            init_code.append('!')
                    # Add sample input
                    if 'code_input' in self[item]:
                        if len(self[item]['code_input']) > 0:
                            init_code.append('! Sample Input')
                            init_code.append('!')
                            for inline in self[item]['code_input']:
                               init_code.append('!  %s'%inline)
                            init_code.append('!')
                    # Add any inline documentation from the input call
                    if 'var_doc' in self[item]:
                        if len(self[item]['var_doc']) > 0:
                            init_code.append('! Input call documentation')
                            init_code.append('!')
                            for inline in self[item]['var_doc']:
                               init_code.append('!  %s'%inline)
                            init_code.append('!')
                            
                    # Add initialization           
                               
                    if 'defvals' in self[item] and 'vars' in self[item]:
                        for var,defval in zip(self[item]['vars'],self[item]['defvals']):
                            init_code.append('   %s = %s'%(var,defval))

                    elif sub == 'rdvmf':
                        init_code.append('  ****************   Special ********************')       

                    else:
                        print("Problem with initialization for TAG=",item)

        #for line in init_code:
        #   print('IC:',line)

        print("Total tags with initialization code:",item_cnt)
        
        with open(name, 'w') as f:
            f.write('\n'.join(init_code))

            

    def create_read_code(self,name):
        # generate the code to read the different series of tags
        read_code = []
        # loop through tagseries
        item_cnt = 0
        for series in self.tagseries:
            
            series_cnt = 0
            read_code.append('')
            read_code.append('!')
            read_code.append('!  Inputs for TAG series %s'%(series))
            read_code.append('!')

            for item in self:
                # ignore comments
                if item.startswith('__comment_'):
                    continue

                # process tags matching the tag series - compose output file
                # contains sample input - commented original read code - converted read code

                # sample input -
                #
                tag = self[item]['tag']
                s = tag[0:1]
                
                if s == series:
                    # start with if statement or elseif statement
                    item_cnt += 1
                    if series_cnt == 0:
                        line = '  if (tag(1:3) .eq. \'%s\') then' % (tag)
                    else:
                        line = '  elseif (tag(1:3) .eq. \'%s\') then' % (tag)

                    series_cnt += 1
                    read_code.append(line)
                    # follow with commented out input examples as documentation
                    line ='!   Sample input '
                    read_code.append(line)
                    if 'code_input' in self[item]:
                        for input in self[item]['code_input']:
                            line ='!   %s' % (input)
                            read_code.append(line)
                    if 'code_comments' in self[item]:
                        for comment in self[item]['code_comments']:
                            read_code.append(comment)
                            
                    # follow with call to load input
                    line = '      call divrd(%s)' % self[item]['args']
                    read_code.append(line)
                    
            read_code.append('  endif')    
            read_code.append('')
            read_code.append('')

        #for line in read_code:
        #    print('RC:',line)

        print("Total tags with initialization code:",item_cnt)        
        
        with open(name, 'w') as f:
            f.write('\n'.join(read_code))

        
    def create_documentation(self,name):
        # combine the input lines to obtain some basic documentation for each option
        print("stub")

        
                        
                    
    def create_output(self):
        # set line lengths to give a more nicely formatted output file that is more legible
        # data for some inputs in OEDGE needs to be within column 72
        # in-line comments can be beyond the end of the line
        tag_width = 4
        comment_width = 50
        end_width = 20
        data_width = 10

        # create output list
        out = []
        for item in self:
            if item.startswith('__comment_'):
                out.append(self[item].strip().ljust(comment_width))

            elif '{' in item:
                out.append('\'%s %s\'' % (item.ljust(tag_width), self[item]['comment'].strip().ljust(comment_width)))
                if isinstance(self[item]['data'], list):
                    for line in self[item]['data']:
                        out.append(' '.join(item.rjust(data_width) for item in map(repr, line)))
                else:
                    out[-1] += ' ' + str(self[item]['data']).strip().rjust(data_width)

            elif 'header_comment' in self[item]:
                if self[item]['comment'] == '':
                    out.append(
                        '\'%s %s\' \'%s\'' % (item.ljust(tag_width), self[item]['comment'].strip(), self[item]['inline_comment'].strip())
                    )
                    for dummy in self[item]['dummies']:
                        out.append(dummy.strip().ljust(comment_width))
                    out.append(
                        '\'%s\' %s %s'
                        % (
                            self[item]['header'].strip().ljust(comment_width),
                            str(len(self[item]['data'])).rjust(data_width),
                            self[item]['header_comment'].rjust(end_width),
                        )
                    )
                    for line in self[item]['data']:
                        out.append(' '.join(item.rjust(data_width) for item in map(repr, line)))
                else:
                    out.append(
                        '\'%s %s\'  %s ' % (item.ljust(tag_width), self[item]['comment'].strip(), self[item]['inline_comment'].strip())
                    )
                    for dummy in self[item]['dummies']:
                        out.append(dummy.strip().ljust(comment_width))
                    out.append(
                        '\'%s\' %s %s'
                        % (
                            self[item]['header'].strip().ljust(comment_width),
                            str(len(self[item]['data'])).rjust(data_width),
                            self[item]['header_comment'].rjust(end_width),
                        )
                    )
                    for line in self[item]['data']:
                        out.append(' '.join(item.rjust(data_width) for item in map(repr, line)))

            elif 'inline_comment' in self[item]:
                if self[item]['comment'] == '':
                    out.append(
                        '\'%s %s\' %s %s'
                        % (
                            item.ljust(tag_width),
                            self[item]['comment'].strip(),
                            ' '.join(item.rjust(data_width) for item in map(repr, tolist(self[item]['data']))),
                            self[item]['inline_comment'],
                        )
                    )
                elif isinstance(self[item]['data'], str):
                    out.append(
                        '\'%s %s\' %s %s'
                        % (
                            item.ljust(tag_width),
                            self[item]['comment'].strip(),
                            ' '.join(item.rjust(data_width) for item in map(repr, tolist(self[item]['data']))),
                            self[item]['inline_comment'],
                        )
                    )
                else:
                    out.append(
                        '\'%s %s\' %s %s'
                        % (
                            item.ljust(tag_width),
                            self[item]['comment'].strip().ljust(comment_width),
                            ' '.join(item.rjust(data_width) for item in map(repr, tolist(self[item]['data']))),
                            self[item]['inline_comment'].rjust(end_width),
                        )
                    )

        return out

    # Use the OMFITObject base class method "saveas" to save using a different file name
    def altsave(self, filename=None):

        if filename is None:
            filename = self.filename

        out = self.create_output()

        with open(filename, 'w') as f:
            f.write('\n'.join(out))

    @dynaSave
    # Use the OMFITObject base class method "saveas" to save using a different file name
    def save(self, **kw):

        out = self.create_output()

        with open(self.filename, 'w') as f:
            f.write('\n'.join(out))


class OMFIToedgeRun(SortedDict):
    pass

    def __init__(self, casename=None, caseloc=None, runloc=None, plot_file_name=None, grid=None, background=None, fc=None, divref=None):
        self.casename = casename
        self.caseloc = caseloc
        self.runloc = runloc
        self.basestr = ''
        self['casename'] = self.casename
        self['plot'] = self.plot_file_name
        self['grid'] = self.casename
        self['casename'] = self.casename

    def __call__(self, *args):
        # Allow the run object to return the string representation of its tree elements
        # by checking case is eval(case()) in the calling routine - can verify object is loaded where expected
        if len(args) == 0:
            return self.basestr
        else:
            tmp = self.basestr
            for arg in args:
                tmp = tmp + "['" + arg + "']"
            return tmp


class OMFIToedgeNC(OMFITnc):
    """OEDGE output NC data file"""

    plot_types_available = ['contour', 'along ring']
    # possible future plots to be added
    # plot_types_available = ['contour','along ring','diagnostic','along wall','along target']

    # background plasma quantities
    # present for cases with a background plasma (basically everything)
    bg_plots = OrderedDict()

    bg_plots['ne'] = {
        'data': ['KNBS'],
        'targ': ['KNDS'],
        'label': 'Density',
        'units': 'm-3',
        'lineaxes': ['P', 'S'],
        'group': False,
        'notebook': False,
    }
    bg_plots['Te'] = {
        'data': ['KTEBS'],
        'targ': ['KTEDS'],
        'label': 'Te',
        'units': 'eV',
        'lineaxes': ['P', 'S'],
        'group': False,
        'notebook': False,
    }
    bg_plots['Ti'] = {
        'data': ['KTIBS'],
        'targ': ['KTIDS'],
        'label': 'Ti',
        'units': 'eV',
        'lineaxes': ['P', 'S'],
        'group': False,
        'notebook': False,
    }
    bg_plots['Vb'] = {
        'data': ['KVHS'],
        'targ': ['KVDS'],
        'label': 'Vpara',
        'units': 'm/s',
        'lineaxes': ['P', 'S'],
        'group': False,
        'notebook': False,
        'scale': '1/QTIM',
    }
    bg_plots['Plasma'] = {
        'data': ['ne', 'Te', 'Ti', 'Vb'],
        'label': 'Background Plasma',
        'lineaxes': ['P', 'S'],
        'group': True,
        'notebook': False,
    }

    bg_plots['E'] = {
        'data': ['KES'],
        'targ': ['KEDS'],
        'label': 'Epara',
        'units': 'V/m',
        'lineaxes': ['P', 'S'],
        'group': False,
        'notebook': False,
    }
    bg_plots['ExB Pol'] = {
        'data': ['E_POL'],
        'targ': None,
        'label': 'Epol',
        'units': 'V/m',
        'lineaxes': ['P', 'S'],
        'group': False,
        'notebook': False,
    }
    bg_plots['ExB Rad'] = {
        'data': ['E_RAD'],
        'targ': None,
        'label': 'Epol',
        'units': 'V/m',
        'lineaxes': ['P', 'S'],
        'group': False,
        'notebook': False,
    }
    bg_plots['Efields'] = {
        'data': ['E', 'ExB Pol', 'ExB Rad'],
        'label': 'Electric Fields',
        'lineaxes': ['P', 'S'],
        'group': True,
        'notebook': False,
    }

    bg_plots['V_pol'] = {
        'data': ['EXB_P'],
        'targ': None,
        'label': 'Vpol',
        'units': 'm/s',
        'lineaxes': ['P', 'S'],
        'group': False,
        'notebook': False,
    }
    bg_plots['V_rad'] = {
        'data': ['EXB_R'],
        'targ': None,
        'label': 'Vrad',
        'units': 'm/s',
        'lineaxes': ['P', 'S'],
        'group': False,
        'notebook': False,
    }
    bg_plots['Drifts'] = {'data': ['V_pol', 'V_rad'], 'label': 'Drifts', 'lineaxes': ['P', 'S'], 'group': True, 'notebook': False}

    # Data present for cases that ran EIRENE to generate hydrogenic quantities
    h_plots = OrderedDict()
    h_plots['H Power'] = {
        'data': ['HPOWLS'],
        'targ': None,
        'label': 'Hydrogen Radiated Power',
        'units': 'W/m3',
        'lineaxes': ['P', 'S'],
        'group': False,
        'notebook': False,
    }
    h_plots['H Line'] = {
        'data': ['HLINES'],
        'targ': None,
        'label': 'Hydrogen Line Radiation',
        'units': 'W/m3',
        'lineaxes': ['P', 'S'],
        'group': False,
        'notebook': False,
    }
    h_plots['H Dalpha'] = {
        'data': ['PINALP'],
        'targ': None,
        'label': 'Hydrogen Dalpha Emission',
        'units': 'ph/m3/s',
        'lineaxes': ['P', 'S'],
        'group': False,
        'notebook': False,
    }
    h_plots['H Ionization'] = {
        'data': ['PINION'],
        'targ': None,
        'label': 'Hydrogen Ionization',
        'units': '/m3/s',
        'lineaxes': ['P', 'S'],
        'group': False,
        'notebook': False,
    }
    h_plots['H Recomb'] = {
        'data': ['PINREC'],
        'targ': None,
        'label': 'Hydrogen Recombination',
        'units': '/m3/s',
        'lineaxes': ['P', 'S'],
        'group': False,
        'notebook': False,
    }
    h_plots['H2 Density'] = {
        'data': ['PINMOL'],
        'targ': None,
        'label': 'Hydrogen Molecule Density',
        'units': '/m3',
        'lineaxes': ['P', 'S'],
        'group': False,
        'notebook': False,
    }
    h_plots['H0 Density'] = {
        'data': ['PINATO'],
        'targ': None,
        'label': 'Hydrogen Atom Density',
        'units': '/m3',
        'lineaxes': ['P', 'S'],
        'group': False,
        'notebook': False,
    }
    h_plots['H Atom Temp'] = {
        'data': ['PINENA'],
        'targ': None,
        'label': 'Hydrogen Atom Temperature',
        'units': 'eV',
        'lineaxes': ['P', 'S'],
        'group': False,
        'notebook': False,
    }
    h_plots['H Mol Temp'] = {
        'data': ['PINENM'],
        'targ': None,
        'label': 'Hydrogen Molecule Temperature',
        'units': 'eV',
        'lineaxes': ['P', 'S'],
        'group': False,
        'notebook': False,
    }
    h_plots['H Ion Energy Loss'] = {
        'data': ['PINQI'],
        'targ': None,
        'label': 'Hydrogen-Ion Energy Loss Term',
        'units': 'W/m3(?)',
        'lineaxes': ['P', 'S'],
        'group': False,
        'notebook': False,
    }
    h_plots['H Elec Energy Loss'] = {
        'data': ['PINQE'],
        'targ': None,
        'label': 'Hydrogen-Electron Energy Loss Term',
        'units': 'W/m3(?)',
        'lineaxes': ['P', 'S'],
        'group': False,
        'notebook': False,
    }
    h_plots['H Quantities'] = {
        'data': ['H0 Density', 'H2 Density', 'H Ionization', 'H Recomb'],
        'label': 'Hydrogen Quantities',
        'lineaxes': ['P', 'S'],
        'group': False,
        'notebook': False,
    }

    # Data present for cases that ran impurities
    imp_plots = OrderedDict()
    imp_plots['Imp Density'] = {
        'data': ['DDLIMS'],
        'targ': None,
        'label': 'Impurity Density',
        'units': '/m3',
        'lineaxes': ['P', 'S'],
        'group': False,
        'notebook': False,
        'scale': 'ABSFAC',
    }
    imp_plots['Imp Temperature'] = {
        'data': ['DDTS'],
        'targ': None,
        'label': 'Impurity Temperature',
        'units': 'eV',
        'lineaxes': ['P', 'S'],
        'group': False,
        'notebook': False,
    }
    imp_plots['Imp Ionization'] = {
        'data': ['TIZS'],
        'targ': None,
        'label': 'Impurity Ionization',
        'units': '/m3/s',
        'lineaxes': ['P', 'S'],
        'group': False,
        'notebook': False,
        'scale': 'ABSFAC',
    }
    imp_plots['Imp Radiated Power'] = {
        'data': ['POWLS'],
        'targ': None,
        'label': 'Impurity Radiated Power',
        'units': 'W/m3',
        'lineaxes': ['P', 'S'],
        'group': False,
        'notebook': False,
        'scale': 'ABSFAC',
    }
    imp_plots['Impurity Quantities'] = {
        'data': ['Imp Density', 'Imp Temperature', 'Imp Ionization', 'Imp Radiated Power'],
        'label': 'Impurity Quantities',
        'lineaxes': ['P', 'S'],
        'group': True,
        'notebook': False,
    }

    # Plots along surfaces (?) separate category or make a different dict entry under appropriate sections?
    surface_plots = OrderedDict()

    # Synthetic diagnostic plots are a separate category since they will require special treatment
    diagnostic_plots = OrderedDict()

    # available plot types - contour/false colour plots, line plots along ring, line plots along surfaces (use different plotting axes)
    plot_types = ['contour-polygon', 'contour-interpolate', 'along ring', 'surface']
    # default_contour = 'contour-polygon'

    all_plots = OrderedDict()

    all_plots.update(bg_plots)
    all_plots.update(h_plots)
    all_plots.update(imp_plots)
    all_plots.update(surface_plots)

    plot_categories = ['background', 'hydrogen', 'impurity', 'surface', 'diagnostic']

    def __init__(self, filename, **kw):
        OMFITnc.__init__(self, filename, **kw)
        self.dynaLoad = True

    @dynaLoad
    def load(self):

        # empty the dictionary in case this is a reload
        self.clear()

        # print('Load netcdf:')

        OMFITnc.load(self)

        #
        # Do some post-processing after loading
        #
        # add anything else it should automatically do with the data on load ...
        # should these be done in the load routine or in the constructor?
        # Presumably the constructor calls load()

        # print('Load self variables:')
        self.load_simulation_data()
        self.plots_available()

    def load_simulation_data(self):
        # This routine labels specific geometry data from the NC file into local variables
        # This data is frequently used in the plots .. since it appears that one case object
        # handles all the plotting then it is worthwhile adding some efficiencies

        # Certain data must be present to make plots possible ... if it is not available for
        # some reason then flag it as an error.

        # load wall and create path object to define inside and outside the vessel
        # This is used in the tricontourf plots to eliminate triangles outside the vessel
        debug = True

        # vessel wall definition
        self.nves = self['NVES']['data']
        self.rves = self['RVES']['data'][: self.nves]
        self.zves = self['ZVES']['data'][: self.nves]
        self.wall = [[self.rves[i], self.zves[i]] for i in range(self.nves)]
        self.wall_path = matplotlib.path.Path(self.wall)

        # load the computational mesh polygon data into a PolyCollection for plotting code results
        # remove zero volume polygons from the list and also exclude when loading simulation data

        self.nvertp = self['NVERTP']['data']
        self.rvertp = self['RVERTP']['data']
        self.zvertp = self['ZVERTP']['data']
        self.korpg = self['KORPG']['data']

        # cell centers
        self.rs = self['RS']['data']
        self.zs = self['ZS']['data']

        # Significant surfaces in the computational mesh
        self.nrs = self['NRS']['data']
        self.nks = self['NKS']['data']
        self.irsep = self['IRSEP']['data']
        # no irsep2 in netcdf
        # self.irsep2 = self['IRSEP2']['data']
        self.irwall = self['IRWALL']['data']
        self.irwall2 = self['IRWALL2']['data']
        self.irtrap = self['IRTRAP']['data']
        self.irtrap2 = self['IRTRAP2']['data']

        # grid connection map data
        self.ikins = self['IKINS']['data']
        self.ikouts = self['IKOUTS']['data']
        self.irins = self['IRINS']['data']
        self.irouts = self['IROUTS']['data']

        # properties of the grid and target
        self.area = self['KAREAS']['data']

        # idds is an index array mapping from (targ,ir) -> (target index)
        # However, the index is from fortran arrays and numbers from 1
        # so when mapped to python arrays you need to subtract 1 which
        # should work if idds is a np array
        self.idds = self['IDDS']['data'].copy() - 1
        self.smax = self['KSMAXS']['data']
        self.pmax = self['KPMAXS']['data']
        self.kss = self['KSS']['data']
        self.kps = self['KPS']['data']

        mesh = []
        count = 0
        debug = True

        for ir in range(self.nrs):
            for ik in range(self.nks[ir]):
                index = self.korpg[ir, ik] - 1
                if self.area[ir, ik] != 0.0:
                    mesh.append(list(zip(self.rvertp[index][0:4], self.zvertp[index][0:4])))
                    count = count + 1

                # check to see that the cell center point is inside the set of vertices for the cell
                # this verifies that the grid is loading correctly
                if debug and self.area[ir, ik] != 0.0:
                    vert = list(zip(self.rvertp[index][0:4], self.zvertp[index][0:4]))
                    cell = matplotlib.path.Path(vert)
                    r = self['RS']['data'][ir, ik]
                    z = self['ZS']['data'][ir, ik]
                    if not cell.contains_point([r, z]):
                        print("grid loading error:")
                        print("vert:", vert)
                        print("test:", ir, ik, r, z, cell.contains_point([r, z]))

        # save the number of non-zero volume cells in the mesh
        # save the mesh polygons as in grid
        self.num_cells = count
        self.grid = mesh

        # The grid is used to create a polygon plot
        # coll = PolyCollection(self.grid,array=data,cmap=pyplot.cm.jet,edgecolors='none')

        #
        # The above created a mesh of polygons from the computational grid
        # Also needed is a mesh composed of the cell center points both with and without the
        # target coordinates. These meshes have two uses:
        # 1) Interpolated triangular meshes which are based on cell center and target data
        # 2) Interpolation of the OEDGE results to generate an ERO plasma file.
        #

        self.num_points = self.num_cells + (self.nrs - self.irsep + 1) * 2

        self.r_cells = np.zeros(self.num_cells)
        self.z_cells = np.zeros(self.num_cells)

        self.r_points = np.zeros(self.num_points)
        self.z_points = np.zeros(self.num_points)

        cp = 0
        cc = 0

        for ir in range(self.nrs):
            for ik in range(self.nks[ir]):

                if ik == 0 and ir >= self.irsep - 1:
                    # add first target point
                    id0 = self.idds[1, ir]
                    self.r_points[cp] = self['RP']['data'][id0]
                    self.z_points[cp] = self['ZP']['data'][id0]
                    cp = cp + 1

                if self.area[ir, ik] != 0.0:
                    self.r_points[cp] = self.rs[ir][ik]
                    self.z_points[cp] = self.zs[ir][ik]
                    cp = cp + 1
                    self.r_cells[cc] = self.rs[ir][ik]
                    self.z_cells[cc] = self.zs[ir][ik]
                    cc = cc + 1

                if ik == self.nks[ir] - 1 and ir >= self.irsep - 1:
                    # add last target point
                    id1 = self.idds[0, ir]
                    self.r_points[cp] = self['RP']['data'][id1]
                    self.z_points[cp] = self['ZP']['data'][id1]
                    cp = cp + 1

        # Impurity scaling factor if available
        # All impurity results are stored scaled to 1 part/m-tor/s entering the system
        # to get absolute values the results are multipled by absfac (absolute scaling factor)
        if 'ABSFAC' in self:
            self.absfac = self['ABSFAC']['data']
        else:
            self.absfac = 1.0

        # Simulation time steps
        # ion time step
        if 'QTIM' in self:
            self.qtim = self['QTIM']['data']
        else:
            self.qtim = 1.0

        # neutral time step
        if 'FSRATE' in self:
            self.fsrate = self['FSRATE']['data']
        else:
            self.fsrate = 1.0

        # print('qtim:',self.qtim,self.fsrate)

        if 'NIZS' in self:
            self.nizs = self['NIZS']['data']
        else:
            self.nizs = None

    def plots_available(self):
        # this routine returns a list of quantities that this object can plot
        # background plasma plots
        self.plots_avail = []
        self.bg_plots_avail = []
        self.h_plots_avail = []
        self.imp_plots_avail = []
        self.surf_plots_avail = []

        for plotkey in list(OMFIToedgeNC.all_plots.keys()):
            # Check to see that all the required datasets are available in the NC file for each class of plots
            avail = True
            for datakey in OMFIToedgeNC.all_plots[plotkey]['data']:
                dataname, targname, label, scalef = self.get_names(plotkey, datakey)
                if dataname not in self:
                    avail = False
                if targname is not None:
                    if targname not in self:
                        avail = False

            # if all elements of a group plot are available then that plot is available
            if avail:
                self.plots_avail.append(plotkey)
                if plotkey in OMFIToedgeNC.bg_plots:
                    self.bg_plots_avail.append(plotkey)
                elif plotkey in OMFIToedgeNC.h_plots:
                    self.h_plots_avail.append(plotkey)
                elif plotkey in OMFIToedgeNC.imp_plots:
                    self.imp_plots_avail.append(plotkey)
                elif plotkey in OMFIToedgeNC.surface_plots:
                    self.surf_plots_avail.append(plotkey)

    def get_plot_categories(self):
        # return a list of the available plot categories
        return self.plot_categories

    def get_plots_available(self, kind=None):
        # return a list of the plots that are available in the current dataset
        # Since the plot list is getting long - allow for splitting by type so that the UI
        # can have more concise drop down lists
        # print("get_plots_available:",kind)
        if kind is None:
            return self.plots_avail
        elif kind == 'background':
            return self.bg_plots_avail
        elif kind == 'hydrogen':
            return self.h_plots_avail
        elif kind == 'impurity':
            return self.imp_plots_avail
        else:
            return []

    def get_plot_types(self):
        # this routine returns a list of the supported plot types
        return self.plot_types

    def get_along_ring_axis_types(self, plot_select):
        # look up the selected plot and return the along ring or linear axis types available
        if plot_select in OMFIToedgeNC.all_plots:
            return OMFIToedgeNC.all_plots[plot_select]['lineaxes']
        # if this fails for some reason then return None ... but it shouldn't
        return None

    def need_ionization_state(self, selection):
        # Check to see if a charge state is needed for the specified plot
        res = False
        if selection in OMFIToedgeNC.imp_plots:
            res = True
        return res

    def need_ring_number(self, selection):
        # Check to see if a ring number is needed for specified plot type
        res = False
        if selection == 'along ring':
            res = True
        return res

    def get_data_2d(self, dataname, targname=None, charge=None, scalef=1.0):
        # load a set of 2D data
        #

        rawdata = self[dataname]['data']
        count = 0
        #
        # Note: targname should always be None for cases where charge is specified
        # charge is a leading index into the data array mapping.

        # if targname is specified then add the target data at the ends of the rings
        if targname is None:
            # leave out target data when plotting a polygon mesh
            data = np.zeros(self.num_cells)
        else:
            # add space to include target values in the mesh - used in triangulation interpolation plots
            data = np.zeros(self.num_points)

        # loop through data arrays
        # Note ... all of this is required because the rows of the array are of inconsistent lengths .. it isn't a rectangular grid
        # in terms of indexing ... core, sol and pfz rings may all have different numbers of cells though nsol_cells=ncore+npfz is usually true
        # along a ring
        count = 0
        for ir in range(self.nrs):
            for ik in range(self.nks[ir]):
                if targname is not None and ik == 0 and ir >= self.irsep - 1:
                    # add first target point
                    id0 = self.idds[1, ir]
                    data[count] = self[targname]['data'][id0]
                    count = count + 1

                if self.area[ir, ik] != 0.0:
                    # include index for charge state if specified - this requires an offset
                    if charge is None:
                        data[count] = rawdata[ir][ik] * scalef
                    else:
                        # Note: the results array indexing goes from -1:nizs for charge states.
                        # -1 maps to 0 in the python array ... charge 0 is in index 1 - so add one to charge
                        # to get correct indexing
                        data[count] = rawdata[charge + 1][ir][ik] * scalef

                    count = count + 1

                if targname is not None and ik == self.nks[ir] - 1 and ir >= self.irsep - 1:
                    # add last target point
                    id1 = self.idds[0, ir]
                    data[count] = self[targname]['data'][id1]
                    count = count + 1

        return data

    def plot(self, plot_select, plot_type, axis_type='S', ring_range=[], charge_range=[], zoom=None):
        # produce a contour plot
        # Either interpolated or polygon

        # print('plot inputs1:',plot_select,plot_type,axis_type,ring_range,charge_range,zoom)

        # Set inputs not needed for specific plots to default values
        # Argument ELIMINATION doesn't propagate ... while argument object MANIPULATION would propagate back to calling routine
        # This is needed since I am not updating the GUI when updating ring or charge ranges ... so these variables keep
        # values from previous plot calls that may not be relevant
        # I could just reload the GUI to avoid this but I am not sure the overhead is worth it
        if not self.need_ring_number(plot_type):
            ring_range = []

        if not self.need_ionization_state(plot_select):
            charge_range = []

        # print('plot inputs2:',plot_select,plot_type,axis_type,ring_range,charge_range,zoom)

        if plot_type == 'contour-polygon' or plot_type == 'contour-interpolate':
            self.plot_contour(plot_select, plot_type, charge_range, zoom)

        elif plot_type == 'along ring':
            if len(ring_range) == 0:
                print("ERROR: Rings not specified correctly for along ring plots")
            else:
                self.plot_along_ring(plot_select, axis_type, ring_range, charge_range, zoom)

        elif plot_type == 'surface':
            print("Plot not yet implemented")
            return

        elif plot_type == 'diagnostic':
            print("Plot not yet implemented")
            return

        else:
            print("ERROR: Specified plot type: `%s` not found." % plot_type)

    def plot_contour(self, plot_select, plot_type, charge_range=[], zoom=None):
        # This produces polygon contour plots matching the
        # computationals mesh of the data included in the plot list

        # If nplots is 1 then the 'data' reference contains a self[] reference otherwise it contains a list of all_plot dictionary entries
        # set up subplots

        #
        # If charge states are specified for an impurity quantity two things happen
        # 1) The figures go into a figure notebook - one/charge state
        # 2) The get_data_2d routine in the specific contour plotting routines will require the charge state to be specified.
        #

        if len(charge_range) == 0:
            # if no charge states then just plot the figure
            # Note: OMFITx.py Figure is intended to put a figure inside an OMFIT GUI element
            # fig = Figure(returnFigure=True)
            fig = figure()
            self.plot_contour_fig(fig, plot_select, plot_type, charge=None, zoom=zoom)
        else:
            # use a figure notebook
            # For lists of charge states
            from omfit_plot import FigureNotebook

            fig_notebook = FigureNotebook(nfig=0, name='Contour plots by charge state')

            # loop through charge states
            for charge in charge_range:
                fig = fig_notebook.add_figure(label='IZ={}'.format(charge))
                self.plot_contour_fig(fig, plot_select, plot_type, charge, zoom=zoom)
                fig.canvas.draw()

    def plot_contour_fig(self, fig, plot_select, plot_type, charge=None, zoom=None):
        # given the figure - add the appropriate contour plots

        plot_list = OMFIToedgeNC.all_plots[plot_select]['data']
        nplts = len(plot_list)

        # get figure layout for number of plots on the page
        nrows, ncols = self.calculate_layout(nplts)

        pyplot.figure(num=fig.number)
        for cnt in range(len(plot_list)):
            ax = pyplot.subplot(nrows, ncols, cnt + 1)
            dataname, targname, label, scalef = self.get_names(plot_select, plot_list[cnt])

            if plot_type == 'contour-polygon':
                # print("plot_contour_polygon_call:",scalef)
                self.plot_contour_polygon(fig, ax, dataname, charge, zoom, scalef=scalef)
                ax.set_title(label)

            elif plot_type == 'contour-interpolate':
                self.plot_contour_interpolate(fig, ax, dataname, targname, charge, zoom, scalef=scalef)
                ax.set_title(label)
        # set the figure title
        fig.suptitle(plot_select)

    def plot_contour_polygon(self, fig, ax, dataname, charge=None, zoom=None, scalef=1.0):
        # polygon plots never use the target data
        # produce polygon contour plot
        if charge is None:
            data = self.get_data_2d(dataname, scalef=scalef)
        else:
            data = self.get_data_2d(dataname, targname=None, charge=charge, scalef=scalef)

        coll = matplotlib.collections.PolyCollection(self.grid, array=data, cmap=pyplot.cm.jet, edgecolors='none')
        ax.add_collection(coll)
        ax.plot(self.rves, self.zves, linewidth=1)
        # ax.plot(wall,linewidth=1)
        ax.autoscale_view()
        ax.set_aspect("equal")
        fig.colorbar(coll, ax=ax)

    def plot_contour_interpolate(self, fig, ax, dataname, targname=None, charge=None, zoom=None, scalef=1.0):
        # produce interpolated contour plot
        # not all 2D data has target data
        if charge is None:
            data = self.get_data_2d(dataname, targname, scalef=scalef)
        else:
            data = self.get_data_2d(dataname, targname, charge, scalef=scalef)

        if targname == None:
            r = self.r_cells
            z = self.z_cells
        else:
            r = self.r_points
            z = self.z_points

        triang = matplotlib.tri.Triangulation(r, z)

        # Mask off unwanted triangles by finding ones with center points outside of wall
        # This could be changed to ones with any corner outside wall but then you might get
        # missing triangle gaps inside the wall
        rmid = r[triang.triangles].mean(axis=1)
        zmid = z[triang.triangles].mean(axis=1)
        points = [[rmid[i], zmid[i]] for i in range(len(rmid))]

        res = self.wall_path.contains_points(points)
        mask = np.where(res == False, 1, 0)
        # print("mask",mask)

        # triangles to be plotted should be set now
        triang.set_mask(mask)

        # ax.set_cmap("jet")
        ax.set_aspect('equal')
        ax.plot(self.rves, self.zves, linewidth=1)
        # contour levels are set to 20 but this could be pulled into a settable parameter
        cs = ax.tricontourf(triang, data, 20, cmap=pyplot.cm.jet)
        fig.colorbar(cs, ax=ax)

    def calculate_layout(self, nplts):
        #
        # Calculate a layout that will hold all the figures
        #
        # Handle special cases - then use generic algorithm
        #
        if nplts == 3:
            nrows = 1
            ncols = 3
        elif nplts == 7:
            nrows = 2
            ncols = 4
        else:
            nrows = np.ceil(np.sqrt(nplts))
            ncols = np.ceil(np.sqrt(nplts))
            # if number of plots will fit in a smaller grid remove a row
            if ncols * (nrows - 1) >= nplts:
                nrows = nrows - 1

        # return the layout
        return nrows, ncols

    def get_names(self, plot_select, plot_name):
        # This returns the names from the NC dictionary to be plotted based on the OMFIToedgeNC.all_plots summary
        # It also returns the label associated with the individual plot element
        dataname = None
        targname = None
        label = None
        scale_desc = None

        # each of the dictionary entries contains a list - this list is either a single reference to
        # the NC datafile or a list of reference to other plots. Plot_name should always be a singular plot if it is present

        if plot_select in OMFIToedgeNC.all_plots:
            if plot_name in self:
                # Need to return plot_name as dataname and get targname from the plot select reference
                dataname = plot_name
                if 'targ' in OMFIToedgeNC.all_plots[plot_select]:
                    if OMFIToedgeNC.all_plots[plot_select]['targ'] is not None:
                        targname = OMFIToedgeNC.all_plots[plot_select]['targ'][0]
                # get label
                if 'label' in OMFIToedgeNC.all_plots[plot_select]:
                    if OMFIToedgeNC.all_plots[plot_select]['label'] is not None:
                        label = OMFIToedgeNC.all_plots[plot_select]['label']
                        if 'units' in OMFIToedgeNC.all_plots[plot_select]:
                            if OMFIToedgeNC.all_plots[plot_select]['units'] is not None:
                                label = label + ' ' + OMFIToedgeNC.all_plots[plot_select]['units']
                # get scaling factor if any
                if 'scale' in OMFIToedgeNC.all_plots[plot_select]:
                    if OMFIToedgeNC.all_plots[plot_select]['scale'] is not None:
                        scale_desc = OMFIToedgeNC.all_plots[plot_select]['scale']

            elif plot_name in OMFIToedgeNC.all_plots:
                # Need to load both dataname and targname from the plot list
                if 'data' in OMFIToedgeNC.all_plots[plot_name]:
                    if OMFIToedgeNC.all_plots[plot_name]['data'] is not None:
                        dataname = OMFIToedgeNC.all_plots[plot_name]['data'][0]

                if 'targ' in OMFIToedgeNC.all_plots[plot_name]:
                    if OMFIToedgeNC.all_plots[plot_name]['targ'] is not None:
                        targname = OMFIToedgeNC.all_plots[plot_name]['targ'][0]
                # get data label
                if 'label' in OMFIToedgeNC.all_plots[plot_name]:
                    if OMFIToedgeNC.all_plots[plot_name]['label'] is not None:
                        label = OMFIToedgeNC.all_plots[plot_name]['label']
                        if 'units' in OMFIToedgeNC.all_plots[plot_name]:
                            if OMFIToedgeNC.all_plots[plot_name]['units'] is not None:
                                label = label + ' ' + OMFIToedgeNC.all_plots[plot_name]['units']
                # get scaling factor if any
                if 'scale' in OMFIToedgeNC.all_plots[plot_name]:
                    if OMFIToedgeNC.all_plots[plot_name]['scale'] is not None:
                        scale_desc = OMFIToedgeNC.all_plots[plot_name]['scale']

        #
        #            Comment out - this code is used both to get plots AND to check whether they exist in the nc file
        #
        #            else:
        #                print("ERROR: plot_name must be either a reference to self or a plot in OMFIToedgeNC.all_plots",plot_select,plot_name)

        # set numerical scaling factor for the volumetric data
        # Scaling factor is set 1.0 by default ... only changed if some other quantity is specified
        if scale_desc == 'ABSFAC':
            scalef = self.absfac
        elif scale_desc == '1/QTIM':
            scalef = 1.0 / self.qtim
        else:
            scalef = 1.0

        return dataname, targname, label, scalef

    def get_data_along_ring(self, ir, dataname=None, targname=None, charge=None, scalef=1.0):
        # loads a set of 2D data
        # get along ring data
        # print("get_data_along_ring:",ir,dataname,targname,charge)

        if dataname is not None:
            if charge is None:
                data = self[dataname]['data'][ir][: self.nks[ir]] * scalef
            else:
                # Note: the results array indexing goes from -1:nizs for charge states.
                # -1 maps to 0 in the python array ... charge 0 is in index 1 - so add one to charge
                # to get correct indexing
                data = self[dataname]['data'][charge + 1][ir][: self.nks[ir]] * scalef

        # print("data:",dataname,data)

        # Insert the target data at the beginning and end when it is available
        # Targt data is only available on sol and pfz rings
        if targname is not None and ir >= self.irsep - 1:
            # Target data does not need the scaling
            # get target indices in the target data arrays
            id0 = self.idds[1, ir]
            id1 = self.idds[0, ir]

            # print('id0,id1:',id0,id1,self.idds[1][ir],self.idds[0][ir])
            # print('idds:',self.idds)

            # print('targ:',targname,self[targname]['data'])
            # print('targ2 id:',targname,self.idds[0][:])
            # print('targ1 id:',targname,self.idds[1][:])
            # print('targ data1:',targname,id0,self[targname]['data'][id0])

            data = np.insert(data, 0, [self[targname]['data'][id0]])

            # print('targ data2:',targname,id1,self[targname]['data'][id1])

            data = np.append(data, [self[targname]['data'][id1]])
            # print('data2:',data)

        # return the data
        return data

    def get_axis(self, ir, axis_type='S', targname=None):
        # This loads the particular axis for the along ring plot (either S or P)
        # It also adds the target elements if targname is not None
        if axis_type == 'P':
            axis = self.kps[ir][: self.nks[ir]]
            label = 'P (m)'
            if targname is not None and ir > self.irsep - 1:
                axis = np.insert(axis, 0, 0.0)
                axis = np.append(axis, self.pmax[ir])
        else:  #  if axis_type == 'S':  # if someone forgets to set axis_type treat it as 'S'
            axis = self.kss[ir][: self.nks[ir]]
            label = 'S (m)'
            if targname is not None and ir >= self.irsep - 1:
                axis = np.insert(axis, 0, 0.0)
                axis = np.append(axis, self.smax[ir])

        return axis, label

    def plot_ring(self, ir, fig, plot_select, axis_type='S', charge_range=[]):
        # plot selected data for specified ring on figure specified

        if ir < 1 or ir > self.nrs:
            print("ERROR: Ring specified is out of range available: %d" % ir)
            return
        else:
            # ring references are in 1 .. nrs from fortran while python arrays index from zero
            ir_ref = ir - 1

        plot_list = OMFIToedgeNC.all_plots[plot_select]['data']
        nplts = len(plot_list)

        # If nplots is 1 then the 'data' reference contains a self[] reference otherwise it contains a list of all_plot dictionary entries

        # set up subplots

        nrows, ncols = self.calculate_layout(nplts)

        pyplot.figure(num=fig.number)
        for cnt in range(len(plot_list)):
            ax = pyplot.subplot(nrows, ncols, cnt + 1)
            if len(charge_range) == 0:
                dataname, targname, label, scalef = self.get_names(plot_select, plot_list[cnt])
                axis, x_label = self.get_axis(ir_ref, axis_type, targname)
                data = self.get_data_along_ring(ir_ref, dataname, targname, scalef=scalef)

                # add a marker for the end points if targname is not None
                if targname is not None:
                    markers_on = [0, len(axis) - 1]
                else:
                    markers_on = []
                # plot the data
                ax.plot(axis, data, linestyle='-', marker='o', markevery=markers_on, linewidth=1)
                ax.set_xlabel(x_label)
                ax.set_ylabel(label)
                ax.set_xlim(0.0, axis[-1])
            else:
                for charge in charge_range:
                    dataname, targname, label, scalef = self.get_names(plot_select, plot_list[cnt])
                    axis, x_label = self.get_axis(ir_ref, axis_type, targname)
                    data = self.get_data_along_ring(ir_ref, dataname, targname, charge, scalef=scalef)
                    # add a marker for the end points if targname is not None
                    if targname is not None:
                        markers_on = [0, len(axis) - 1]
                    else:
                        markers_on = []
                    # plot the data
                    ax.plot(axis, data, label='IZ={}'.format(charge), linestyle='-', marker='o', markevery=markers_on, linewidth=1)
                    ax.set_xlabel(x_label)
                    ax.set_ylabel(label)
                    ax.set_xlim(0.0, axis[-1])
                ax.legend()

    def plot_along_ring(self, plot_select, axis_type, ring_range, charge_range=[], zoom=None):
        # Routine plots the quantities along ring ... number of plots/page is controlled by number
        # of items in the plot list.

        # Single figures
        # Single ring/single quantity
        # single ring/multiple quantities
        #
        # Figure notebook spanning ring range
        # multi rings/single quantity
        # multi rings/multiple quantities
        #
        # If more than one ring is specified then a figure notebook is used
        #
        # If more than one charge - try to put them on the same axes
        #
        # plot

        nrings = len(ring_range)
        if nrings <= 0:
            # error condition
            print("ERROR: No rings specified for the plot")
            return

        minring = 1
        maxring = self.nrs

        # print('rings:',nrings,ring_range)

        res = [ring_range[i] < minring or ring_range[i] > maxring for i in range(len(ring_range))]

        # print('res:',res)

        out_of_bounds = np.where(res == True)
        # print('out_of_bounds:',len(out_of_bounds),)
        # print('test:',ring_range[out_of_bounds])

        # Check to see if at least some rings are in range
        if len(out_of_bounds) >= nrings:
            print("ERROR: No valid rings numbers specified: min=%d max=%d range=%s" % (minring, maxring, str(ring_range)))
            return

        plot_list = OMFIToedgeNC.all_plots[plot_select]['data']
        nplts = len(plot_list)

        # if nplts is 1 then the list contains a self reference ... if greater than 1 then it contains a list of plot references
        # Do everything in a figure notebook since there can be 1+ figures generated

        # print('rings:',nrings,ring_range,plot_select,axis_type)

        # use figure notebook since it supports an unknown number of plots
        from omfit_plot import FigureNotebook

        fig_notebook = FigureNotebook(nfig=0, name='Along Ring Plots')

        for ir in ring_range:
            fig = fig_notebook.add_figure(label='IR={}'.format(ir))
            self.plot_ring(ir, fig, plot_select, axis_type, charge_range)
            fig.canvas.draw()
