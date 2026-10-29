// Manuscript layout: A4, modest margins, title block with lettered affiliations, compact floats.
#let article(
  title: none,
  subtitle: none,
  authors: none,
  affiliations: none,
  date: none,
  abstract: none,
  abstract-title: none,
  cols: 1,
  margin: (x: 2.2cm, y: 2.4cm),
  paper: "a4",
  lang: "en",
  region: "US",
  font: "libertinus serif",
  fontsize: 10.5pt,
  title-size: 1.5em,
  subtitle-size: 1.25em,
  heading-family: "libertinus serif",
  heading-weight: "bold",
  heading-style: "normal",
  heading-color: black,
  heading-line-height: 0.65em,
  sectionnumbering: none,
  pagenumbering: "1",
  toc: false,
  toc_title: none,
  toc_depth: none,
  toc_indent: 1.5em,
  doc,
) = {
  set page(
    paper: paper,
    margin: margin,
    numbering: pagenumbering,
  )
  set par(justify: true, leading: 0.62em, spacing: 0.9em)
  set text(lang: lang, region: region, font: font, size: fontsize)
  set heading(numbering: sectionnumbering)
  show heading: set block(above: 1.5em, below: 0.8em, sticky: true)
  set bibliography(title: none)
  show heading.where(level: 1): set text(size: 1.2em)
  show heading.where(level: 2): set text(size: 1.05em)
  show figure.where(kind: image): it => place(auto, float: true, clearance: 1.6em, it)
  show figure.caption: set text(size: 0.88em)
  show figure.caption: set par(justify: true)
  show table: set text(size: 0.82em)
  show link: set text(hyphenate: false)
  show raw: set text(size: 0.9em)

  if title != none {
    align(center)[
      #block(above: 0pt, below: 1.4em)[
        #set par(leading: heading-line-height, justify: false)
        #text(weight: "bold", size: title-size)[#title]
        #if subtitle != none {
          parbreak()
          text(weight: "bold", size: subtitle-size)[#subtitle]
        }
      ]
    ]
  }

  if authors != none {
    align(center)[
      #block(below: 0.9em)[
        #set par(justify: false)
        #authors.map(a => [#a.name#super[#a.affiliations.map(n => numbering("a", n)).join(",")#if a.corresponding [,\*]]]).join([, ], last: [ and ])
      ]
      #if affiliations != none {
        block(below: 0.6em)[
          #set par(justify: false, leading: 0.5em)
          #set text(size: 0.9em)
          #for aff in affiliations [
            #super[#numbering("a", aff.number)] #aff.text \
          ]
          #for a in authors.filter(a => a.corresponding) [
            #super[\*]Corresponding author: #a.name, #a.email
          ]
        ]
      }
    ]
  }

  if date != none {
    align(center)[#block(inset: 1em)[#date]]
  }

  if abstract != none {
    block(inset: (x: 1em, y: 1em))[
      #text(weight: "semibold")[#abstract-title] #h(1em) #abstract
    ]
  }

  if toc {
    let title = if toc_title == none { auto } else { toc_title }
    block(above: 0em, below: 2em)[
      #outline(title: toc_title, depth: toc_depth, indent: toc_indent);
    ]
  }

  if cols == 1 {
    doc
  } else {
    columns(cols, doc)
  }
}

#set table(
  inset: 5pt,
  stroke: none
)
