%%% Modification of BibTeX style file dcu.bst
%%% ... by urlbst, version 0.5 (marked with "% urlbst")
%%% See <http://www.astro.gla.ac.uk/users/norman/distrib/urlbst/>
%%% Added webpage entry type, and url and lastchecked fields.
%%% Original headers follow...

% BibTeX standard bibliography style `dcu' (one of the harvard family)
	% version 0.99a for BibTeX versions 0.99a or later, LaTeX version 2.09.
	% Copyright (C) 1991, all rights reserved.
	% Copying of this file is authorized only if either
	% (1) you make absolutely no changes to your copy, including name, or
	% (2) if you do make changes, you name it something other than
	% btxbst.doc, plain.bst, unsrt.bst, alpha.bst, abbrv.bst, agsm.bst,
        % dcu.bst or kluwer.bst.
	% This restriction helps ensure that all standard styles are identical.
	% The file harvard.tex has the documentation for this style.

% ACKNOWLEDGEMENT:
%   This document is a modified version of alpha.bst to which it owes much of
%   its functionality.

% AUTHOR
%   Peter Williams, Key Centre for Design Quality, Sydney University
%   e-mail: peterw@archsci.arch.su.oz.au

ENTRY
  { address
    author
    booktitle
    chapter
    edition
    editor
    howpublished
    institution
    journal
    key
    month
    note
    number
    organization
    pages
    publisher
    school
    series
    title
    type
    volume
    year
    eprint % urlbst
    doi % urlbst
    url % urlbst
    lastchecked % urlbst
  }
  { field.used }
  { extra.label sort.label list.year }

INTEGERS { output.state before.all mid.sentence after.sentence after.block }

STRINGS { urlintro eprinturl eprintprefix doiprefix doiurl } % urlbst...
INTEGERS { hrefform addeprints adddoiresolver }
% Following constants may be adjusted by hand, if desired
FUNCTION {init.config.constants}
{
  "Available from: " 'urlintro := % prefix before URL
  "http://arxiv.org/abs/" 'eprinturl := % prefix to make URL from eprint ref
  "arXiv:" 'eprintprefix := % text prefix printed before eprint ref
  "http://dx.doi.org/" 'doiurl := % prefix to make URL from DOI
  "doi:" 'doiprefix := % text prefix printed before DOI ref
  #0 'addeprints := % 0=no eprints; 1=include eprints
  #0 'adddoiresolver := % 0=no DOI resolver; 1=include it
  #0 'hrefform := % 0=no crossrefs; 1=hypertex xrefs; 2=hyperref refs
}
INTEGERS { 
  bracket.state
  outside.brackets
  open.brackets
  within.brackets
  close.brackets
}
% ...urlbst to here
FUNCTION {init.state.consts}
{ #0 'outside.brackets := % urlbst
  #1 'open.brackets :=
  #2 'within.brackets :=
  #3 'close.brackets :=

  #0 'before.all :=
  #1 'mid.sentence :=
  #2 'after.sentence :=
  #3 'after.block :=
}

STRINGS { s t f }

% urlbst
FUNCTION {output.nonnull.original}
{ 's :=
  output.state mid.sentence =
    { ", " * write$ }
    { output.state after.block =
	{ add.period$ write$
	  newline$
	  "\newblock " write$
	}
	{ output.state before.all =
	    'write$
	    { add.period$ " " * write$ }
	  if$
	}
      if$
      mid.sentence 'output.state :=
    }
  if$
  s
}

% urlbst...
FUNCTION {output.nonnull}
{ % Save the thing we've been asked to output
  's :=
  % If the bracket-state is close.brackets, then add a close-bracket to
  % what is currently at the top of the stack, and set bracket.state
  % to outside.brackets
  bracket.state close.brackets =
    { "]" *
      outside.brackets 'bracket.state :=
    }
    'skip$
  if$
  bracket.state outside.brackets =
    { % We're outside all brackets -- this is the normal situation.
      % Write out what's currently at the top of the stack, using the
      % original output.nonnull function.
      s
      output.nonnull.original
    }
    { % Still in brackets.  Add open-bracket or (continuation) comma, add the
      % new text (in s) to the top of the stack, and move to the close-brackets
      % state, ready for next time (unless inbrackets resets it).  If we come
      % into this branch, then output.state is carefully undisturbed.
      bracket.state open.brackets =
        { " [" * }
        { ", " * } % bracket.state will be within.brackets
      if$ 
      s * 
      close.brackets 'bracket.state :=
    }
  if$
}

% Call this function just before adding something which should be presented in 
% brackets.  bracket.state is handled specially within output.nonnull.
FUNCTION {inbrackets}
{ bracket.state close.brackets =
    { within.brackets 'bracket.state := } % reset the state: not open nor closed
    { open.brackets 'bracket.state := }
  if$
}

FUNCTION {format.lastchecked}
{ lastchecked empty$
    { "" }
    { inbrackets "cited " lastchecked * }
  if$
}
% ...urlbst to here

FUNCTION {output}
{ duplicate$ empty$
    'pop$
    'output.nonnull
  if$
}

FUNCTION {output.check}
{ 't :=
  duplicate$ empty$
    { pop$ "empty " t * " in " * cite$ * warning$ }
    'output.nonnull
  if$
}

FUNCTION {item.check}
{ 't :=
  empty$
    { "empty " t * " in " * cite$ * warning$ }
    { skip$ }
  if$
}

FUNCTION {fin.entry}
{ 
  bracket.state close.brackets = % urlbst
    { "]" * }
    'skip$
  if$
   add.period$
  write$
  newline$
}

FUNCTION {new.block}
{ output.state before.all =
    'skip$
    { after.block 'output.state := }
  if$
}

FUNCTION {not}
{   { #0 }
    { #1 }
  if$
}

FUNCTION {and}
{   'skip$
    { pop$ #0 }
  if$
}

FUNCTION {or}
{   { pop$ #1 }
    'skip$
  if$
}

FUNCTION {field.or.null}
{ duplicate$ empty$
    { pop$ "" }
    'skip$
  if$
}

FUNCTION {emphasize}
{ duplicate$ empty$
    { pop$ "" }
    { "{\em " swap$ * "}" * }
  if$
}

FUNCTION {embolden}
{ duplicate$ empty$
    { pop$ "" }
    { "{\textbf " swap$ * "}" * }
  if$
}

FUNCTION {quote}
{ duplicate$ empty$
    { pop$ "" }
    { "`" swap$ * "'" * }
  if$
}

INTEGERS { nameptr namesleft numnames }

FUNCTION {format.names}
{ 's :=
  'f :=
  #1 'nameptr :=
  s num.names$ 'numnames :=
  numnames 'namesleft :=
    { namesleft #0 > }
    { s nameptr f format.name$ 't :=
      nameptr #1 >
	{ namesleft #1 >
	    { ", " * t * }
	    { t "others" =
		{ " et~al." * }
		{ " and " * t * }
	      if$
	    }
	  if$
	}
	't
      if$
      nameptr #1 + 'nameptr :=
      namesleft #1 - 'namesleft :=
    }
  while$
}

FUNCTION {format.authors}
{ author empty$
    { "" }
    { "{vv~}{ll}{, jj}{, f.}" author format.names }
  if$
}

FUNCTION {format.editors}
{ editor empty$
    { "" }
    { "{vv~}{ll}{, jj}{, f.}" editor format.names
      editor num.names$ #1 >
	{ " (eds)" * }
	{ " (ed.)" * }
      if$
    }
  if$
}

FUNCTION {format.editors.reverse}
{ editor empty$
    { "" }
    { "{f.~}{vv~}{ll}{, jj}" editor format.names
      editor num.names$ #1 >
	{ " (eds)" * }
	{ " (ed.)" * }
      if$
    }
  if$
}

FUNCTION {format.title}
{ title empty$
    { "" }
    { title "t" change.case$ }
  if$
}

FUNCTION {n.dashify}
{ 't :=
  ""
    { t empty$ not }
    { t #1 #1 substring$ "-" =
	{ t #1 #2 substring$ "--" = not
	    { "--" *
	      t #2 global.max$ substring$ 't :=
	    }
	    {   { t #1 #1 substring$ "-" = }
		{ "-" *
		  t #2 global.max$ substring$ 't :=
		}
	      while$
	    }
	  if$
	}
	{ t #1 #1 substring$ *
	  t #2 global.max$ substring$ 't :=
	}
      if$
    }
  while$
}

FUNCTION {format.btitle}
{ title emphasize
}

FUNCTION {tie.or.space.connect}
{ duplicate$ text.length$ #3 <
    { "~" }
    { " " }
  if$
  swap$ * *
}

FUNCTION {either.or.check}
{ empty$
    'pop$
    { "can't use both " swap$ * " fields in " * cite$ * warning$ }
  if$
}

FUNCTION {format.bvolume}
{ volume empty$
    { "" }
    { "Vol." volume tie.or.space.connect
      series empty$
	'skip$
	{ " of " * series emphasize * }
      if$
      "volume and number" number either.or.check
    }
  if$
}

FUNCTION {format.number.series}
{ volume empty$
    { number empty$
	{ series field.or.null }
	{ output.state mid.sentence =
	    { "number" }
	    { "Number" }
	  if$
	  number tie.or.space.connect
	  series empty$
	    { "there's a number but no series in " cite$ * warning$ }
	    { " in " * series emphasize * }
	  if$
	}
      if$
    }
    { "" }
  if$
}

FUNCTION {format.edition}
{ edition empty$
    { "" }
    { output.state mid.sentence =
	{ edition "l" change.case$ " edn" * }
	{ edition "t" change.case$ " edn" * }
      if$
    }
  if$
}

INTEGERS { multiresult }

FUNCTION {multi.page.check}
{ 't :=
  #0 'multiresult :=
    { multiresult not
      t empty$ not
      and
    }
    { t #1 #1 substring$
      duplicate$ "-" =
      swap$ duplicate$ "," =
      swap$ "+" =
      or or
	{ #1 'multiresult := }
	{ t #2 global.max$ substring$ 't := }
      if$
    }
  while$
  multiresult
}

FUNCTION {format.pages}
{ pages empty$
    { "" }
    { pages multi.page.check
	{ "pp.~" pages n.dashify * }
	{ "p.~" pages * }
      if$
    }
  if$
}

FUNCTION {format.vol.num.pages}
{ volume embolden field.or.null
  number empty$
    'skip$
    { "(" number * ")" * *
      volume empty$
	{ "there's a number but no volume in " cite$ * warning$ }
	'skip$
      if$
    }
  if$
  pages empty$
    'skip$
    { duplicate$ empty$
	{ pop$ format.pages }
	{ ":~" * pages n.dashify * }
      if$
    }
  if$
}

FUNCTION {format.chapter.pages}
{ chapter empty$
    'format.pages
    { type empty$
	{ "chapter" }
	{ type "l" change.case$ }
      if$
      chapter tie.or.space.connect
      pages empty$
	'skip$
	{ ", " * format.pages * }
      if$
    }
  if$
}

FUNCTION {format.in.ed.booktitle}
{ booktitle empty$
    { "" }
    { editor empty$
	{ booktitle emphasize }
	{ "{\em in} " format.editors.reverse * ", " * booktitle emphasize * }
      if$
    }
  if$
}

FUNCTION {empty.misc.check}
{ author empty$ title empty$ howpublished empty$
  month empty$ year empty$ note empty$
  and and and and and
  key empty$ not and
    { "all relevant fields are empty in " cite$ * warning$ }
    'skip$
  if$
}

FUNCTION {format.thesis.type}
{ type empty$
    'skip$
    { pop$
      type "t" change.case$
    }
  if$
}

FUNCTION {format.tr.number}
{ type empty$
    { "Technical Report" }
    'type
  if$
  number empty$
    { "t" change.case$ }
    { number tie.or.space.connect }
  if$
}

FUNCTION {format.article.crossref}
{ key empty$
    { journal empty$
	{ "need key or journal for " cite$ * " to crossref " * crossref *
	  warning$
	  ""
	}
	{ "in {\em " journal * "\/}" * " \cite{" * crossref * "}" *
        }
      if$
    }
    { " {\em in} \citeasnoun{" crossref * "}" * }
  if$
}

FUNCTION {format.book.crossref}
{ volume empty$
    { "empty volume in " cite$ * "'s crossref of " * crossref * warning$
      "in "
    }
    { "Vol." volume tie.or.space.connect
      " of " *
    }
  if$
  editor empty$
  editor field.or.null author field.or.null =
  or
    { key empty$
	{ series empty$
	    { "need editor, key, or series for " cite$ * " to crossref " *
	      crossref * warning$
	      "" *
	    }
	    { "{\em " * series * "\/}" * " \cite{" * crossref * "}" *}
	  if$
	}
	{ " \citeasnoun{" * crossref * "}" * }
      if$
    }
    { " \citeasnoun{" * crossref * "}" * }
  if$
}

FUNCTION {format.incoll.inproc.crossref}
{ editor empty$
  editor field.or.null author field.or.null =
  or
    {
      key empty$
	{ booktitle empty$
	    { "need editor, key, or booktitle for " cite$ * " to crossref " *
	      crossref * warning$
	      ""
	    }
	    { "in {\em " booktitle * "\/}" * " \cite{" * crossref * "}" *}
	  if$
	}
	{ " {\em in} \citeasnoun{" crossref * "}" * }
      if$
    }
    {  " {\em in} \citeasnoun{" crossref * "}" * }
  if$
}

INTEGERS { len }

FUNCTION {chop.word}
{ 's :=
  'len :=
  s #1 len substring$ =
    { s len #1 + global.max$ substring$ }
    's
  if$
}

FUNCTION {format.lab.names.abbr}
{ 's :=
  s num.names$ 'numnames :=
  numnames #1 >
    { numnames #2 >
	{ s #1 "{vv~}{ll}" format.name$ " {\textit and others}" * }
	{ s #2 "{ff }{vv }{ll}{ jj}" format.name$ "others" =
            { s #1 "{vv~}{ll}" format.name$ " {\textit and others}" * }
	    { s #1 "{vv~}{ll}" format.name$ " and " *
              s #2 "{vv~}{ll}" format.name$ * 
            }
          if$
        }
      if$
    }
    { s #1 "{vv~}{ll}" format.name$ }
  if$
}

FUNCTION {format.lab.names.full}
{ 's :=
  #1 'nameptr :=
  s num.names$ 'numnames :=
  numnames 'namesleft :=
    { namesleft #0 > }
    { s nameptr "{vv~}{ll}" format.name$ 't :=
      nameptr #1 >
	{ namesleft #1 >
	    { ", " * t * }
	    { t "others" =
		{ " et~al." * }
		{ " and " * t * }
	      if$
	    }
	  if$
	}
	't
      if$
      nameptr #1 + 'nameptr :=
      namesleft #1 - 'namesleft :=
    }
  while$
}

INTEGERS { author.field editor.field organization.field title.field key.field }

FUNCTION {init.field.constants}
{ #0 'author.field :=
  #1 'editor.field :=
  #2 'organization.field :=
  #3 'title.field :=
  #4 'key.field :=
}

FUNCTION {make.list.label}
{ author.field field.used =
    { format.authors }
    { editor.field field.used =
        { format.editors }
        { organization.field field.used =
            { "The " #4 organization chop.word #3 text.prefix$ }
            { title.field field.used =
                { format.btitle }
                { key.field field.used =
                    { key #3 text.prefix$ }
                    { "Internal error :001 on " cite$ * " label" * warning$ }
                  if$
                }
              if$
            }
          if$
        }
      if$
    }
  if$
}

FUNCTION {make.full.label}
{ author.field field.used =
    { author format.lab.names.full }
    { editor.field field.used =
        { editor format.lab.names.full }
        { organization.field field.used =
            { "The " #4 organization chop.word #3 text.prefix$ }
            { title.field field.used =
                { format.btitle }
                { key.field field.used =
                    { key #3 text.prefix$ }
                    { "Internal error :001 on " cite$ * " label" * warning$ }
                  if$
                }
              if$
            }
          if$
        }
      if$
    }
  if$
}

FUNCTION {make.abbr.label}
{ author.field field.used =
    { author format.lab.names.abbr }
    { editor.field field.used =
        { editor format.lab.names.abbr }
        { organization.field field.used =
            { "The " #4 organization chop.word #3 text.prefix$ }
            { title.field field.used =
                { format.btitle }
                { key.field field.used =
                    { key #3 text.prefix$ }
                    { "Internal error :001 on " cite$ * " label" * warning$ }
                  if$
                }
              if$
            }
          if$
        }
      if$
    }
  if$
}

FUNCTION {output.bibitem}
{ outside.brackets 'bracket.state := % urlbst
   newline$
  "\harvarditem[" write$
  make.abbr.label write$
  "]{" write$
  make.full.label write$
  "}{" write$
  list.year write$
  "}{" write$
  cite$ write$
  "}" write$
  newline$
  ""
  before.all 'output.state :=
}

FUNCTION {list.label.output}
{ make.list.label " " * write$
}

FUNCTION {format.date}
{ year empty$
    { month empty$
	{ "" }
	{ "there's a month but no year in " cite$ * warning$
	  month
	}
      if$
    }
    { month empty$
	'year
	{ month " " * year * }
      if$
    }
  if$
}
% urlbst...
% Functions for making hypertext links.
% In all cases, the stack has (link-text href-url)
%
% make 'null' specials
FUNCTION {make.href.null}
{
  pop$
}
% make hypertex specials
FUNCTION {make.href.hypertex}
{ 
  "\special {html:<a href=" quote$ *
  swap$ * quote$ * "> }" * swap$ *
  "\special {html:</a>}" *
}
% make hyperref specials
FUNCTION {make.href.hyperref}
{ 
  "\href {" swap$ * "} {" * swap$ * "}" *
}
FUNCTION {make.href}
{ hrefform #2 =
    'make.href.hyperref      % hrefform = 2
    { hrefform #1 =
        'make.href.hypertex  % hrefform = 1
        'make.href.null      % hrefform = 0 (or anything else)
      if$
    }
  if$
}

FUNCTION {format.url}
{ url empty$
    { "" }
      { hrefform #1 =
          { % special case -- add HyperTeX specials
            urlintro "\url{" url * "}" * url make.href.hypertex * }
          { urlintro "\url{" * url * "}" * }
       if$
     }
  if$
}

FUNCTION {format.eprint}
{ eprint empty$
    { "" }
    { eprintprefix eprint * eprinturl eprint * make.href }
  if$
}

FUNCTION {format.doi}
{ doi empty$
    { "" }
    { doiprefix doi * doiurl doi * make.href }
  if$
}

% Output a URL.  We can't use the more normal idiom (something like
% `format.url output'), because the `inbrackets' within
% format.lastchecked applies to everything between calls to `output',
% so that `format.url format.lastchecked * output' ends up with both
% the URL and the lastchecked in brackets.
FUNCTION {output.url}
{ url empty$
    'skip$ 
    { new.block 
      format.url output
      format.lastchecked output 
    }
  if$
}

FUNCTION {output.web.refs}
{
  new.block
  output.url
  addeprints eprint empty$ not and
    { format.eprint output.nonnull }
    'skip$
  if$
  adddoiresolver doi empty$ not and
    { format.doi output.nonnull }
    'skip$
  if$
%  addeprints
%    { eprint empty$
%        'skip$
%        { format.eprint output.nonnull }
%      if$
%    }
%    'skip$
%  if$
}

% Webpage entry type.
% Title and url fields required;
% author, note, year, month, and lastchecked fields optional
% See references 
%   ISO 690-2 http://www.nlc-bnc.ca/iso/tc46sc9/standard/690-2e.htm
%   http://www.classroom.net/classroom/CitingNetResources.html
%   http://neal.ctstateu.edu/history/cite.html
%   http://www.cas.usf.edu/english/walker/mla.html
% for citation formats for web pages.
FUNCTION {webpage}
{ output.bibitem
  author empty$
    { editor empty$
        'skip$  % author and editor both optional
        { format.editors output.nonnull }
      if$
    }
    { editor empty$
        { format.authors output.nonnull }
        { "can't use both author and editor fields in " cite$ * warning$ }
      if$
    }
  if$
%  author empty$
%    'skip$
%    { format.authors output.nonnull }
%  if$
  new.block
  format.title "title" output.check
  inbrackets "online" output
  new.block
  year empty$
    'skip$
    { format.date "year" output.check }
  if$
  lastchecked empty$
    'skip$
    { format.lastchecked output }
  if$
  new.block
  format.url "url" output.check
  new.block
  note output
  fin.entry
}
% ...urlbst to here


FUNCTION {article}
{ output.bibitem
  list.label.output
  " (" list.year * ")" * output.nonnull
  new.block
  author "author" item.check
  title.field field.used =
    { skip$ }
    { format.title "title" output.check }
  if$
  crossref missing$
    { journal emphasize "journal" duplicate$ item.check
      pages empty$
        {
          output
        }
        {
          " " *
          format.vol.num.pages * output
        }
      if$
    }
    { format.article.crossref output.nonnull
      format.pages output
    }
  if$
  new.block
  note output
  output.web.refs  % urlbst
  fin.entry
}

FUNCTION {book}
{ output.bibitem
  list.label.output
  " (" list.year * ")" * output.nonnull
  new.block
  author empty$
    { editor "author and editor" item.check }
    { crossref missing$
	{ "author and editor" editor either.or.check }
	'skip$
      if$
    }
  if$
  title.field field.used =
    { skip$ }
    { format.btitle "title" output.check }
  if$
  crossref missing$
    { format.bvolume output
      format.number.series output
      format.edition output
      publisher "publisher" output.check
      address output
    }
    { format.book.crossref output.nonnull
      format.edition output
    }
  if$
  new.block
  note output
  output.web.refs  % urlbst
  fin.entry
}

FUNCTION {booklet}
{ output.bibitem
  list.label.output
  " (" list.year * ")" * output.nonnull
  new.block
  title.field field.used =
    { skip$ }
    { format.title "title" output.check }
  if$
  howpublished output
  address output
  new.block
  note output
  output.web.refs  % urlbst
  fin.entry
}

FUNCTION {inbook}
{ output.bibitem
  list.label.output
  " (" list.year * ")" * output.nonnull
  new.block
  author empty$
    { editor "author and editor" item.check }
    { crossref missing$
	{ "author and editor" editor either.or.check }
	'skip$
      if$
    }
  if$
  title.field field.used =
    { skip$ }
    { format.btitle "title" output.check }
  if$
  crossref missing$
    { format.bvolume output
      format.number.series output
      format.edition output
      publisher "publisher" output.check
      address output
    }
    { format.book.crossref output.nonnull
      format.edition output
    }
  if$
  format.chapter.pages "chapter and pages" output.check
  new.block
  note output
  output.web.refs  % urlbst
  fin.entry
}

FUNCTION {incollection}
{ output.bibitem
  list.label.output
  " (" list.year * ")" * output.nonnull
  new.block
  title.field field.used =
    { skip$ }
    { format.title "title" output.check }
  if$
  author "author" item.check
  crossref missing$
    { format.in.ed.booktitle "booktitle" output.check
      format.edition output
      format.bvolume output
      format.number.series output
      publisher "publisher" output.check
      address output
    }
    { format.incoll.inproc.crossref output.nonnull
    }
  if$
  format.chapter.pages output
  new.block
  note output
  output.web.refs  % urlbst
  fin.entry
}

FUNCTION {inproceedings}
{ output.bibitem
  list.label.output
  " (" list.year * ")" * output.nonnull
  new.block
  title.field field.used =
    { skip$ }
    { format.title "title" output.check }
  if$
  author "author" item.check
  crossref missing$
    { format.in.ed.booktitle "booktitle" output.check
      format.bvolume output
      format.number.series output
      address empty$
	{ organization output
	  publisher output
	}
	{ organization output
	  publisher output
	  address output.nonnull
	}
      if$
    }
    { format.incoll.inproc.crossref output.nonnull
    }
  if$
  format.pages output
  new.block
  note output
  output.web.refs  % urlbst
  fin.entry
}

FUNCTION {conference} { inproceedings }

FUNCTION {manual}
{ output.bibitem
  list.label.output
  " (" list.year * ")" * output.nonnull
  new.block
  title.field field.used =
    { skip$ }
    { format.btitle "title" output.check }
  if$
  format.edition output
  author empty$
    { organization empty$
	{ address output
	}
	'skip$
      if$
    }
    { organization output
      address output
    }
  if$
  new.block
  note output
  output.web.refs  % urlbst
  fin.entry
}

FUNCTION {mastersthesis}
{ output.bibitem
  list.label.output
  " (" list.year * ")" * output.nonnull
  new.block
  author "author" item.check
  title.field field.used =
    { skip$ }
    { format.title emphasize "title" output.check }
  if$
  "Master's thesis" format.thesis.type output.nonnull
  school "school" output.check
  address output
  new.block
  note output
  output.web.refs  % urlbst
  fin.entry
}

FUNCTION {misc}
{ output.bibitem
  list.label.output
  " (" list.year * ")" * output.nonnull
  new.block
  title.field field.used =
    { skip$ }
    { format.title output }
  if$
  howpublished output
  new.block
  note output
  output.web.refs  % urlbst
  fin.entry
  empty.misc.check
}

FUNCTION {phdthesis}
{ output.bibitem
  list.label.output
  " (" list.year * ")" * output.nonnull
  new.block
  author "author" item.check
  title.field field.used =
    { skip$ }
    { format.btitle "title" output.check }
  if$
  "PhD thesis" format.thesis.type output.nonnull
  school "school" output.check
  address output
  new.block
  note output
  output.web.refs  % urlbst
  fin.entry
}

FUNCTION {proceedings}
{ output.bibitem
  list.label.output
  " (" list.year * ")" * output.nonnull
  new.block
  title.field field.used =
    { skip$ }
    { format.btitle "title" output.check }
  if$
  format.bvolume output
  format.number.series output
  address empty$
    { editor empty$
	{ skip$ }
	{ organization output
	}
      if$
      publisher output
    }
    { editor empty$
	'skip$
	{ organization output }
      if$
      publisher output
      address output.nonnull
    }
  if$
  new.block
  note output
  output.web.refs  % urlbst
  fin.entry
}

FUNCTION {techreport}
{ output.bibitem
  list.label.output
  " (" list.year * ")" * output.nonnull
  new.block
  author "author" item.check
  title.field field.used =
    { skip$ }
    { format.title "title" output.check }
  if$
  format.tr.number emphasize output.nonnull
  institution "institution" output.check
  address output
  new.block
  note output
  output.web.refs  % urlbst
  fin.entry
}

FUNCTION {unpublished}
{ output.bibitem
  list.label.output
  " (" list.year * ")" * output.nonnull
  new.block
  author "author" item.check
  title.field field.used =
    { skip$ }
    { format.title "title" output.check }
  if$
  note "note" output.check
  output.web.refs  % urlbst
  fin.entry
}

FUNCTION {default.type} { misc }

MACRO {jan} {"January"}

MACRO {feb} {"February"}

MACRO {mar} {"March"}

MACRO {apr} {"April"}

MACRO {may} {"May"}

MACRO {jun} {"June"}

MACRO {jul} {"July"}

MACRO {aug} {"August"}

MACRO {sep} {"September"}

MACRO {oct} {"October"}

MACRO {nov} {"November"}

MACRO {dec} {"December"}

MACRO {acmcs} {"ACM Computing Surveys"}

MACRO {acta} {"Acta Informatica"}

MACRO {cacm} {"Communications of the ACM"}

MACRO {ibmjrd} {"IBM Journal of Research and Development"}

MACRO {ibmsj} {"IBM Systems Journal"}

MACRO {ieeese} {"IEEE Transactions on Software Engineering"}

MACRO {ieeetc} {"IEEE Transactions on Computers"}

MACRO {ieeetcad}
 {"IEEE Transactions on Computer-Aided Design of Integrated Circuits"}

MACRO {ipl} {"Information Processing Letters"}

MACRO {jacm} {"Journal of the ACM"}

MACRO {jcss} {"Journal of Computer and System Sciences"}

MACRO {scp} {"Science of Computer Programming"}

MACRO {sicomp} {"SIAM Journal on Computing"}

MACRO {tocs} {"ACM Transactions on Computer Systems"}

MACRO {tods} {"ACM Transactions on Database Systems"}

MACRO {tog} {"ACM Transactions on Graphics"}

MACRO {toms} {"ACM Transactions on Mathematical Software"}

MACRO {toois} {"ACM Transactions on Office Information Systems"}

MACRO {toplas} {"ACM Transactions on Programming Languages and Systems"}

MACRO {tcs} {"Theoretical Computer Science"}

READ

EXECUTE {init.field.constants}

FUNCTION {sortify}
{ purify$
  "l" change.case$
}

FUNCTION {author.key.label}
{ author empty$
    { key empty$
	{ title.field 'field.used := }
	{ key.field 'field.used := }
      if$
    }
    { author.field 'field.used := }
  if$
}

FUNCTION {author.editor.key.label}
{ author empty$
    { editor empty$
	{ key empty$
	    { title.field 'field.used := }
	    { key.field 'field.used := }
	  if$
  	}
	{ editor.field 'field.used := }
      if$
    }
    { author.field 'field.used := }
  if$
}

FUNCTION {author.key.organization.label}
{ author empty$
    { key empty$
	{ organization empty$
	    { title.field 'field.used := }
	    { organization.field 'field.used := }
	  if$
	}
	{ key.field 'field.used := }
      if$
    }
    { author.field 'field.used := }
  if$
}

FUNCTION {editor.key.organization.label}
{ editor empty$
    { key empty$
	{ organization empty$
	    { title.field 'field.used := }
	    { organization.field 'field.used := }
	  if$
	}
	{ key.field 'field.used := }
      if$
    }
    { editor.field 'field.used := }
  if$
}

FUNCTION {sort.format.title}
{ 't :=
  "A " #2
    "An " #3
      "The " #4 t chop.word
    chop.word
  chop.word
  sortify
  #1 global.max$ substring$
}

FUNCTION {calc.label}
{ type$ "book" =
  type$ "inbook" =
  or
    'author.editor.key.label
    { type$ "proceedings" =
	'editor.key.organization.label
	{ type$ "manual" =
	    'author.key.organization.label
	    'author.key.label
	  if$
	}
      if$
    }
  if$
  make.abbr.label
  title.field field.used =
    { sort.format.title }
    { sortify }
  if$
  year field.or.null purify$ #-1 #4 substring$ sortify
  *
  'sort.label :=
}

FUNCTION {first.presort}
{ calc.label
  sort.label
  title.field field.used =
    { skip$ }
    { "    "
      *
      make.list.label sortify
      *
      "    "
      *
      title field.or.null
      sort.format.title
      *
    }
  if$
  #1 entry.max$ substring$
  'sort.key$ :=
}

ITERATE {first.presort}

SORT

STRINGS { last.sort.label next.extra }

INTEGERS { last.extra.num }

FUNCTION {initialize.last.extra.num}
{ #0 int.to.chr$ 'last.sort.label :=
  "" 'next.extra :=
  #0 'last.extra.num :=
}

FUNCTION {forward.pass}
{ last.sort.label sort.label =
    { last.extra.num #1 + 'last.extra.num :=
      last.extra.num int.to.chr$ 'extra.label :=
    }
    { "a" chr.to.int$ 'last.extra.num :=
      "" 'extra.label :=
      sort.label 'last.sort.label :=
    }
  if$
}

FUNCTION {reverse.pass}
{ next.extra "b" =
    { "a" 'extra.label := }
    'skip$
  if$
  year empty$
    { "n.d." extra.label * 'list.year := }
    { year extra.label * 'list.year := }
  if$
  extra.label 'next.extra :=
}

EXECUTE {initialize.last.extra.num}

ITERATE {forward.pass}

REVERSE {reverse.pass}

FUNCTION {second.presort}
{ make.list.label
  title.field field.used =
    { sort.format.title }
    { sortify }
  if$
  "    "
  *
  list.year field.or.null sortify
  *
  "    "
  *
  title.field field.used =
    { skip$ }
    { title field.or.null
      sort.format.title
      *
    }
  if$
  #1 entry.max$ substring$
  'sort.key$ :=
}

ITERATE {second.presort}

SORT

FUNCTION {begin.bib}
{ preamble$ empty$
    'skip$
    { preamble$ write$ newline$ }
  if$
  "\begin{thebibliography}{xx}" write$ newline$
}

EXECUTE {begin.bib}

EXECUTE {init.config.constants}
EXECUTE {init.state.consts}

ITERATE {call.type$}

FUNCTION {end.bib}
{ newline$
  "\end{thebibliography}" write$ newline$
}

EXECUTE {end.bib}
