# Shape Morphing
Morphs shapes on a 2D canvas in a similar way to Adobe After Effects, Flash or GSAP.

## Options
Options are specified using URL parameters.

Name|Function|
--- | ---
vertices | The number of vertices one of the shapes has.
vertices2 | The number of vertices the other shape has.
render | The kind of output to produce. Defaults to animation. Use the keyword `overlay` to produce a static image with the animation frames rendered on top of one another.
speed | The animation speed in pixels/frame.
max_rotation | The maximum amount of rotation which can be applied to transform one shape into the other as a number of degrees between 0 and 180. Larger values produce more coherent transformations where points move together as a whole shape and less like a collection of independent vertex movements. However, smaller values may make you feel less dizzy. Default: 90
fill | The fill colour as a [CSS colour value](https://developer.mozilla.org/en-US/docs/Web/CSS/Guides/Colors/Color_values). When specifying colour parameters `x` can be used as time varying variable which gets substituted with values between 0 and 1 which represent the beginning and end of the animation respectively. Because certain characters have special meanings within URLs you need to use % to escape them. Use %25 to write %, %2B to write +, and %23 to write #.
blend | The [compositing mode](https://developer.mozilla.org/en-US/docs/Web/API/CanvasRenderingContext2D/globalCompositeOperation).
fill2 | A second fill colour. When two fill colours are specified the shape is filled with a colour gradient.
gradient | The method used to calculate the colour gradient. See the table below for permissible values.
stroke | The colour used to draw the shape's outline.
stroke2 | Used to apply a colour gradient to the shape's outline.
stroke_gradient | The method used to calculate the colour gradient for the shape's outline. 2 and 5 are the only algorithms supported.
line_width | The width of the stroke used to paint the shape's outline.
line_width2 | The width of the stroke in the last frame of the animation if it's different from at the start of the animation. If the shape grows or shrinks during the animation then the line width will vary automatically in a way which is proportional to the change to the shape's overall size and it's not necessary to explicitly specify a second line width for this purpose.
dash | A comma separated list of numbers which specify a pattern of dashes and gaps used when drawing the outline.
dash_offset | Whereabouts to begin within the dash pattern.
dash2, dash_offset2 | The dashes can grow, shrink, or shift position during the course of the animation.
start | The shape's outline can be partially drawn. This parameter determines where to begin drawing the shape's outline. It's a value between 0 and 1 which is measured a fraction of the total perimeter distance.
end | A value between 0 and 1 which determines where to stop drawing the shape's outline.
start2, end2 | The section of the outline that's drawn can vary over time can vary.
shadow, shadow2 | The length of the shape's shadow.
shadow_blur | The amount of blurring to be applied to the shadow as a proportion of the shadow's length. This determines how sharply defined the shadow is.
shadow_colour | The colour the shadow is drawn in. This is usually black but other colours can produce interesting artistic effects when combined with the overlay rendering mode.
sun_x, sun_y | The positioning of the light source used to generate the shadow, as proportions of the width and height of the window.
blur | The amount of blur effect to apply to the whole shape, in pixels.
blur2 | The amount of blur to apply in the last frame of the animation (if different from the start of the animation).

Each feature which can animated has its own easing function associated with it which determines if change happens linearly or if it happens with some acceleration at the beginning of the motion and/or deceleration toward the end. A selection of different easing functions have been constructed but the choice of which easing function gets applied to each one of the various visual characteristics cannot be altered from the URL.

### Gradient Types

Type | Description
--- | ---
1 (default) | The gradient begins from somewhere along an edge and runs in a perpendicular direction (if possible) to a vertex opposite the edge.
2 | The gradient runs along one of the edges.
3 | The gradient runs directly from one vertex to another vertex.
4 | The gradient is a [conic gradient](https://developer.mozilla.org/en-US/docs/Web/API/CanvasRenderingContext2D/createConicGradient) with different colours emanating out from a vertex.
5 | The gradient is a conic gradient centred on the centre of the shape. An additional parameter called `spin` can be used to animate the gradient by specifying an amount of rotation in degrees.

## Bugs
* If the shape is concave then the movement paths of two vertices can sometimes cross one another and produce a weird looking result.
* The dashes move around even when I don't intend them to as the shape's perimeter changes.

## Things Not Supported (Yet)
* Shapes with curved parts.
* Shapes with holes, e.g. an annulus (doughnut).
* Shapes with self intersecting edges. E.g. a pentagram.
* Multiple source and destination shapes which would be automatically paired up to minimize the total amount of movement. If there aren't equal numbers of source and destination shapes then some shapes would need to broken into pieces.
* No explicit support/assistance yet with chaining multiple successive animations together.
* A shear transformation isn't currently considered as a possible cause for a shape change. Presently, translations, rotations and scaling of the source shape are analyzed as being possible mappings between the source shape and the destination shape with any residual differences being interpreted as translations of individual vertices.
* Automatic conversion of text to shapes.
* Radial gradients.
* Images used as fill patterns.
* Animation by altering the positions of colours within a gradient rather than changing the colours themselves.
* More alternate rendering modes. For example, a storyboard view containing multiple static images, or a movement path trace for each vertex.
* More easing functions.
* Fractal edges.
