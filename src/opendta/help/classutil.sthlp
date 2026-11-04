{smcl}
{title:Title}

{p2colset 5 22 24 2}{...}
{p2col :{cmd:classutil} {hline 2}}Class programming: .class files, objects and the classutil utility{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}
{cmd:.}{it:name} {cmd:=} {cmd:.}{it:classname}{cmd:.new} [{it:arguments}]{break}
{cmd:.}{it:name}{cmd:.}{it:member} {cmd:=} {it:rvalue}{break}
{cmd:.}{it:name}{cmd:.ref =} {cmd:.}{it:other}{cmd:.ref}{break}
{cmd:.}{it:name}{cmd:.}{it:program} [{it:arguments}]{break}
{cmd:class exit} [{it:rvalue}]{break}
{cmd:classutil drop} {it:instance} [{it:instance} ...]{break}
{cmd:classutil describe} {it:object} [{cmd:,} {cmd:recurse} {cmd:newok}]{break}
{cmd:classutil dir} [{it:pattern}] [{cmd:,} {cmd:all} {cmd:detail}]{break}
{cmd:classutil cdir} [{it:pattern}]{break}
{cmd:classutil which} {it:classname} [{cmd:,} {cmd:all}]

{p 4 4 2}
{cmd:cutil} is a synonym for {cmd:classutil}.

{title:Description}

{pstd}
A class is defined in the file {it:classname}{cmd:.class}, found along the adopath like an
ado-file. The file has a {cmd:version} line, a {cmd:class} {it:classname} {cmd:{c -(}} ... {cmd:{c )-}}
block with the member variables, optionally followed by {cmd:, inherit(}{it:classnames}{cmd:)},
and {cmd:program .}{it:name} ... {cmd:end} blocks with the member programs.
{p_end}

{pstd}
Member variables are {cmd:double}, {cmd:string}, {cmd:array} or an instance of another class,
declared as {it:type name} [{cmd:=} {it:rvalue}] or {it:name} {cmd:=} {it:rvalue}. Variables after
{cmd:classwide:} are shared by all instances; those after {cmd:instancespecific:} (the default)
belong to each instance. Initial values are computed once, when the definition is read.
{p_end}

{pstd}
An {it:rvalue} is {cmd:"}{it:string}{cmd:"}, a number, an expression (in parentheses when it
starts with a quote or a period), a reference {cmd:.}{it:id}[{cmd:.}{it:id} ...] (a copy of the
object, or what a member program returns with {cmd:class exit}), {cmd:{c -(}{c )-}} or
{cmd:{c -(}}{it:el}{cmd:,} {it:el} ...{cmd:{c )-}} for an array. Assigning a value of another type
to an existing object is error {cmd:r(109)}.
{p_end}

{pstd}
Inside a member program, {cmd:.x} means the member {cmd:x} of the instance the program runs on
when it has one, and the global {cmd:.x} otherwise; {cmd:.Global.x} and {cmd:.Local.x} choose
one of them. {cmd:.Super.}{it:prog} (or {cmd:.Super(}{it:classname}{cmd:).}{it:prog}) runs the
program of a parent class. A {cmd:.new} program runs after an instance is created,
{cmd:.oncopy} after a copy (with {cmd:.oncopy_src}, the key of the source), {cmd:.destructor}
before an instance is destroyed, and {cmd:.macroexpand} gives the printable form of an instance.
{p_end}

{pstd}
The built-ins are {cmd:.new}, {cmd:.copy}, {cmd:.ref}, {cmd:.ref_n}, {cmd:.objtype}, {cmd:.isa},
{cmd:.classname}, {cmd:.isofclass} {it:classname}, {cmd:.objkey}, {cmd:.uname}, {cmd:.arrnels},
{cmd:.arrindexof} {cmd:"}{it:string}{cmd:"}, {cmd:.classmv}, {cmd:.instancemv}, {cmd:.dynamicmv},
{cmd:.superclass}, and the modifiers {cmd:.Declare}, {cmd:.Arrdropel}, {cmd:.Arrdropall},
{cmd:.Arrpop} and {cmd:.Arrpush}.
{p_end}

{pstd}
In any command, {cmd:`.}{it:id}{cmd:'} is replaced by the printable form of the object: the
string, the number in {cmd:%18.0g}, or nothing for arrays and instances (unless
{cmd:.macroexpand} is defined). An error in the reference gives nothing.
{p_end}

{pstd}
{cmd:classutil drop} drops top-level instances; {cmd:describe} shows an object or a class;
{cmd:dir} lists the top-level instances (with {cmd:all}, the loaded classes too) and stores
{cmd:r(list)}; {cmd:cdir} lists the {cmd:.class} files on the adopath; {cmd:which} shows the
file of a class and its {cmd:*!} lines and stores {cmd:r(fn)}. {cmd:discard} and
{cmd:clear all} drop every instance and class definition.
{p_end}

{title:Example}

{phang2}{cmd:. .coord = .coordinate.new}{p_end}
{phang2}{cmd:. .coord.x = 1}{p_end}
{phang2}{cmd:. display "`.coord.x'"}{p_end}
{phang2}{cmd:. classutil describe .coord}{p_end}

{title:Note}

{pstd}
The layout of {cmd:classutil describe} and {cmd:classutil dir} is OpenDTA's own, and so is the
form of the keys and unames.
{p_end}
