# Colour plan across Figures 1–5

Written 29 Aug 2026, after "they are all screaming so much red and blue".

## The diagnosis

Red and blue are not the problem. **Saturation is.** Almost every coloured mark
in the current set sits at full strength, so nothing has priority — the eye has
no way to tell the finding from the background. The AASLD heatmap reads calm
not because it uses different hues (it uses the same ones) but because it uses a
*continuous ramp with white at zero*, so most cells are pale and only the
extremes carry weight.

A second, related problem: in several panels colour repeats what position
already shows, so it is spending the loudest channel on nothing.

## Two rules

**1. Colour carries what position does not.**
A bar above zero is already "increased"; a box above the dashed line is already
"up". Colouring those by direction is decoration. Colour earns its place where
it is the *only* encoding — heatmap cells, scatter points, pathway-map nodes.

**2. Saturation is a signal, not a default.**
Full-strength red or blue should mark the few things that carry the finding.
Everything else is a tint, a neutral, or a point on a ramp. Where colour must
carry a quantity, use a continuous diverging ramp with white at zero.

## Audit and proposed changes

| panel | colour now | verdict | change |
|---|---|---|---|
| 1a, 1c volcano | class hues + grey | colour is the only class encoding | keep |
| 1b, 1d class counts | full red / full blue | direction already shown by which side of zero the bar sits | **tint** the fills |
| 2a–c heatmaps | continuous blue–white–red | already a ramp, already calm | keep |
| 2d boxplots | grey / blue | muted, groups not direction | keep |
| 2e saturation bars | purple sequential ramp | ordered variable, ramp is correct | keep |
| 2f index dots | grey / blue | calm | keep |
| 3a boxes, 3b bars | red on grey, sparse | red marks the few enriched sets — this is rule 2 working | keep |
| 3c triglyceride marks | full red / full blue, 264 × 4 marks | colour is the only direction encoding, but at this density it dominates the page | **tint**, keep shape |
| 4a, 4b arm bars | full red / full blue | direction already shown by side; and it competes with panel c | **neutral grey**, so c owns the colour in that row |
| 4c integration | binary → now continuous ramp | fixed 29 Aug | done |
| 4 pathview | tints of the Fig 3 hues | already tinted for label legibility | keep |
| 4e Cyp boxplots | grey / blue | calm | keep |

## What this buys

Figure 4 becomes a hierarchy rather than a wall: a and b are quiet rankings,
c is the coloured payoff that shows where the arms agree, the pathview is a
tinted map, and the boxplots are grey and blue. One loud thing per figure.

Figure 3 keeps its red asterisked sets loud and quiets the 1,056 triangles
behind them.

## Not changing

The hues themselves. Red = increased and blue = decreased is established across
the manuscript and in the reviewers' heads; NC_RED and NC_BLUE stay exactly as
they are. Everything above is about *strength*, not hue.

## Implementation

`nc_style.R` already has `nc_tint()`. Add `NC_RED_SOFT` / `NC_BLUE_SOFT` beside
the existing tints and point the panels above at them. One place to change if
the balance turns out wrong.
