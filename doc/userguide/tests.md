# test markdown
## hyperlinks
[one with a title](http://fsf.org "click here for a good time!"). Unclear how to enforce new window.

## Code environment
Either use fenced style (tildes)

~~~~~~~
if (a > 3) {
  moveShip(5 * gravity, DOWN);
}
~~~~~~~

or indented style (4 whitespaces)

    if (a > 3) {
      moveShip(5 * gravity, DOWN);
    }

Both works with pandoc and wordpress. Also see [pandoc verbatim code](http://pandoc.org/README.html#verbatim-code-blocks "pandoc verbatim code").

## Equations
(@gleichung1) $$a=b*c$$
As (@gleichung1) shows, blabla.

## Bibtex, cite
Hindenlang [@hindenlang2014mesh]. Only works with pandoc!

[bibshow file=https://www.flexi-project.org/wp-content/uploads/2016/07/userguide-1.bib]

Hindenlang [bibcite key=hindenlang2014mesh], Gassner [bibcite key=gassner2011disp]


## section references
## Figures, caption
![This is the caption\label{mylabel}](https://www.flexi-project.org/wp-content/uploads/2016/01/M7_ROE_N7M10_q_0000060p2000000.jpg)

See figure \ref{mylabel}.

## tables
## unnumbered section headings
  just add

    {-}

 after the heading

## Code blocks for various languages

``` {.C}

int a = 32;
int a = 32;

```
