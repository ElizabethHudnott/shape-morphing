import Canvas from './src/canvas.js';
import ColourMorph from './src/colour.js';
import Ease from './src/easing.js';
import {Geometry} from './src/math.js';
import {Line, Shape, RandomShape} from './src/shapes.js';

function alignLeftOrTop(values) {
	const numValues = values.length;
	let min = Infinity;
	for (let i = 0; i < numValues; i++) {
		min = Math.min(min, values[i]);
	}
	for (let i = 0; i < numValues; i++) {
		values[i] -= min;
	}
}

function alignRightOrBottom(values, alignmentValue) {
	const numValues = values.length;
	let max = -Infinity;
	for (let i = 0; i < numValues; i++) {
		max = Math.max(max, values[i]);
	}
	const offset = alignmentValue - max;
	for (let i = 0; i < numValues; i++) {
		values[i] += offset;
	}
}

function centreValues(values, centre) {
	const numValues = values.length;
	let min = Infinity, max = -Infinity;
	for (let i = 0; i < numValues; i++) {
		const value = values[i];
		if (value < min) {
			min = value;
		}
		if (value > max) {
			max = value;
		}
	}
	const offset = centre - min - 0.5 * (max - min);

	for (let i = 0; i < numValues; i++) {
		values[i] += offset;
	}
}

function layoutPolygons(pointsGenerator, numPoints, numPoints2 = numPoints) {
	let width1, height1, width2, height2;
	let pointsX1, pointsY1, pointsX2, pointsY2;
	const width = canvas.width;
	const height = canvas.height;
	const swap = Math.random() < 0.5;
	if (swap) {
		[numPoints, numPoints2] = [numPoints2, numPoints];
	}
	const r = Math.random();
	if (r < 0.5) {
		// Case 1: Top left and bottom right
		// Case 2: Bottom left and top right
		width1 = (0.5 * (3 + Math.random()) / 3) * width;
		height1 = Math.min(width1 / (Math.random() + 1), 2/3 * height);
		width2 = (0.5 * (3 + Math.random()) / 3) * width;
		height2 = Math.min(width2 / (Math.random() + 1), 2/3 * height);
		[pointsX1, pointsY1] = pointsGenerator(numPoints, 0.5 * width1, 0.5 * height1);
		[pointsX2, pointsY2] = pointsGenerator(numPoints2, 0.5 * width2, 0.5 * height2);
		alignLeftOrTop(pointsX1);
		alignRightOrBottom(pointsX2, width);
		if (r < 0.25) {
			alignLeftOrTop(pointsY1);
			alignRightOrBottom(pointsY2, height);
		} else {
			alignLeftOrTop(pointsY2);
			alignRightOrBottom(pointsY1, height);
		}
	} else {
		// Case 3: Entire left and middle right
		// Case 4: Middle left and entire right
		height1 = height;
		width1 = Math.min(height1 / (Math.random() + 1), 0.5 * width);
		width2 = (0.5 * (3 + Math.random()) / 3) * width;
		height2 = Math.min(width2 / (Math.random() + 1), 0.5 * height);
		[pointsX1, pointsY1] = pointsGenerator(numPoints, 0.5 * width1, 0.5 * height1);
		[pointsX2, pointsY2] = pointsGenerator(numPoints2, 0.5 * width2, 0.5 * height2);
		centreValues(pointsY1, 0.5 * height);
		centreValues(pointsY2, 0.5 * height);
		if (r < 0.75) {
			alignLeftOrTop(pointsX1);
			alignRightOrBottom(pointsX2, width);
		} else {
			alignRightOrBottom(pointsX1, width);
			alignLeftOrTop(pointsX2);
		}
	}
	const shape1 = new Shape(pointsX1, pointsY1);
	const shape2 = new Shape(pointsX2, pointsY2);
	return swap ? [shape2, shape1] : [shape1, shape2];
}

function moveLength(polygon1, polygon2) {
	const numPoints = polygon2.numPoints;
	const inverseScale = 1 / polygon2.scaleFactor;

	const translateX = polygon2.centreX - polygon1.centreX;
	const translateY = polygon2.centreY - polygon1.centreY;
	const absRotation = Math.abs(polygon2.rotation);
	let maxTranslation = 0, maxArcLength = 0;
	if (absRotation > 0) {

		for (let i = 0; i < numPoints; i++) {
			const radius = polygon2.radii[i];
			const scaledRadius = radius * inverseScale;
			const radialGrowth = scaledRadius - radius;

			maxTranslation = Math.max(
				maxTranslation,
				Math.abs(translateX + polygon2.offsetsX[i]),
				Math.abs(translateY + polygon2.offsetsY[i])
			);

			// Arc length of an Archimedean spiral.
			const radialGrowthPerRad = radialGrowth / absRotation;
			let squareRoot = Math.sqrt(
				scaledRadius * scaledRadius + radialGrowthPerRad * radialGrowthPerRad
			);
			const firstTerm =  0.5 * squareRoot * scaledRadius / radialGrowthPerRad;
			const secondTerm = 0.5 * radialGrowthPerRad * Math.log(squareRoot + scaledRadius);
			squareRoot = Math.sqrt(radius * radius + radialGrowthPerRad * radialGrowthPerRad);
			const thirdTerm = 0.5 * squareRoot * radius / radialGrowthPerRad;
			const fourthTerm = 0.5 * radialGrowthPerRad * Math.log(squareRoot + radius);
			const arcLength = (firstTerm + secondTerm) - (thirdTerm + fourthTerm);
			maxArcLength = Math.max(maxArcLength, arcLength);
		}

	} else {

		for (let i = 0; i < numPoints; i++) {
			const radius = polygon2.radii[i];
			const scaledRadius = radius * inverseScale;
			const radialGrowth = scaledRadius - radius;
			const angle = polygon2.angles[i];
			const growthX = radialGrowth * Math.cos(angle);
			const growthY = radialGrowth * Math.sin(angle);

			maxTranslation = Math.max(
				maxTranslation,
				Math.abs(translateX + growthX + polygon2.offsetsX[i]),
				Math.abs(translateY + growthY + polygon2.offsetsY[i])
			);
		}

	}
	return Math.max(maxTranslation, maxArcLength);
}

class StrokeStyle {
	constructor() {
		this.width = undefined;
		this.dash = undefined;
		this.dashOffset = 0;
		this.start = 0;
		this.end = 1;
	}
}

class StrokeMorph {

	constructor(closed, colourMorph) {
		this.closed = closed;
		this.colourMorph = colourMorph;
		this.startStyle = new StrokeStyle();
		this.endStyle = new StrokeStyle();

		this.colour = undefined;
		this.width = undefined;
		// TODO scale the dash and dash offset as the perimeter length changes.
		// TODO scale the dash offset by the number of frames?
		this.dash = undefined;
		this.dashOffset = undefined;
		this.start = 0;
		this.end = 1;

		this.widthEase = Ease.linear;
		this.dashOffsetEase = Ease.linear;
		this.startEase = Ease.linear;
		this.endEase = Ease.linear;
	}

	interpolate(morph, interpolation, context) {
		let t;
		if (this.colourMorph) {
			this.colour = this.colourMorph.interpolate(morph, interpolation, context);
		} else {
			this.colour = undefined;
		}

		const startWidth = this.startStyle.width;
		if (startWidth === undefined) {
			this.width = undefined;
		} else {
			const endWidth = this.endStyle.width;
			if (endWidth === undefined) {
				this.width = startWidth;
			} else {
				t = this.widthEase(interpolation);
				this.width =
					(startWidth * (1 - t) + this.endStyle.width * t) /
					(0.5 * (morph.scaleX + morph.scaleY));
			}
		}

		const startDash = this.startStyle.dash;
		if (startDash === undefined) {
			this.dash = undefined;
		} else {
			const numDashTerms = startDash.length;
			const endDash = this.endStyle.dash;
			const dash = new Array(numDashTerms);
			for (let i = 0; i < numDashTerms; i++) {
				dash[i] = startDash[i] * (1 - interpolation) + endDash[i] * interpolation;
			}
			this.dash = dash;
			const startOffset = this.startStyle.dashOffset;
			const endOffset = this.endStyle.dashOffset;
			t = this.dashOffsetEase(interpolation);
			this.dashOffset = startOffset * (1 - t) + endOffset * t;
		}

		const startStart = this.startStyle.start;
		const endStart = this.endStyle.start;
		t = this.startEase(interpolation);
		this.start =  startStart * (1 - t) + endStart * t;

		const startEnd = this.startStyle.end;
		const endEnd = this.endStyle.end;
		t = this.endEase(interpolation);
		this.end =  startEnd * (1 - t) + endEnd * t;
	}

	setColour(colourMorph) {
		this.colourMorph = colourMorph;
	}

	setWidth(start, end = start) {
		this.startStyle.width = start;
		this.endStyle.width = end;
	}

	/**
	 * @param {number[]} start A dash pattern or undefined to leave as is. Use [] to force a
	 * solid line.
	 */
	setDash(start, end = start) {
		if (start === undefined) {
			this.startStyle.dash = undefined;
			return;
		}

		let startLength = start.length;
		let endLength = end.length;
		if (startLength !== endLength) {
			if (startLength & 1) {
				start = start.concat(start);
				startLength <<= 1;
			}
			if (endLength & 1) {
				end = end.concat(end);
				endLength <<= 1;
			} else if (endLength === 0) {
				end = new Array(startLength);
				for (let i = 0; i < startLength; i += 2) {
					end[i] = start[i] + start[i + 1];
					end[i + 1] = 0;
				}
			}
			if (startLength === 0) {
				start = new Array(endLength);
				for (let i = 0; i < endLength; i += 2) {
					start[i] = end[i] + end[i + 1];
					start[i + 1] = 0;
				}
			}
		}

		const commonLength = Geometry.lcm(startLength, endLength);
		for (let i = startLength; i < commonLength; i++) {
			start[i] = start[i % startLength];
		}
		for (let i = endLength; i < commonLength; i++) {
			end[i] = end[i % endLength];
		}
		this.startStyle.dash = start;
		this.endStyle.dash = end;
	}

	setDashOffset(start, end = start) {
		this.startStyle.dashOffset = start;
		this.endStyle.dashOffset = end;
	}

	setStart(startFraction, endFraction = startFraction) {
		this.startStyle.start = startFraction;
		this.endStyle.start = endFraction;
	}

	setEnd(startFraction, endFraction = startFraction) {
		this.startStyle.end = startFraction;
		this.endStyle.end = endFraction;
	}

}

class ShadowMorph {

	constructor(startLength, endLength, blur = 0.25, colourMorph = ColourMorph.Constant.BLACK) {
		this.startLength = startLength;
		this.endLength = endLength;
		this.lengthEase = Ease.quadInOut;
		this.blurFraction = blur;
		this.colourMorph = colourMorph;

		this.offsetX = 0;
		this.offsetY = 0;
		this.blur = 0;
		this.colour = 'black';
	}

	interpolate(morph, interpolation, lightSourceX, lightSourceY) {
		const t = this.lengthEase(interpolation);
		let length = this.startLength * (1 - t) + this.endLength * t;
		length += (morph.strokeMorph?.width || 0) * 0.5;
		const numPoints = morph.numPoints;
		const translateX = morph.translateX, translateY = morph.translateY;
		let sumOffsetX = 0, sumOffsetY = 0;
		for (let i = 0; i < numPoints; i++) {
			let dx = lightSourceX - translateX - morph.pointsX[i];
			let dy = lightSourceY - translateY - morph.pointsY[i];
			if (dx !== 0 || dy !== 0) {
				dx /= morph.scaleX;
				dy /= morph.scaleY;
				const hypotenuse = Math.hypot(dx, dy);
				const offsetX = -dx / hypotenuse * length;
				const offsetY = -dy / hypotenuse * length;
				sumOffsetX += offsetX;
				sumOffsetY += offsetY;
			}
		}
		this.offsetX = sumOffsetX / numPoints;
		this.offsetY = sumOffsetY / numPoints;
		this.blur = this.blurFraction * length;
		this.colour = this.colourMorph.interpolate(morph, interpolation);
	}

}

class StringArtMorph {

	constructor(
		numPoints1, offset1, increment, closed, numPoints2 = numPoints1, offset2 = offset1
	) {
		this.numPoints1 = numPoints1;
		this.offset1 = offset1;
		this.increment = increment % numPoints1;
		this.closed = closed;
		this.numPoints2 = numPoints2;
		this.offset2 = offset2;

		this.startWidth = undefined;
		this.endWidth = undefined;
		this.widthEase = Ease.linear;

		this.colourMorphs = undefined;
		this.colourLengths = undefined;

		this.lines = [];
		this.width = undefined;
	}

	setWidth(start, end = start) {
		this.startWidth = start;
		this.endWidth = end;
	}

	setColours(colourMorphs, lengths) {
		this.colourMorphs = colourMorphs;
		this.colourLengths = lengths;
	}

	interpolate(morph, interpolation, context) {
		const polygon1 = morph.polygon1;
		const polygon2 = morph.polygon2;
		const perimeter1 = Geometry.perimeter(polygon1.deltaX, polygon1.deltaY, this.closed);
		const perimeter2 = Geometry.perimeter(polygon2.deltaX, polygon2.deltaY, this.closed);
		const perimeter = Geometry.perimeter(
			morph.pointsX.map(x => x * morph.scaleX),
			morph.pointsY.map(y => y * morph.scaleY),
			this.closed
		);
		let t = Math.min(Math.max(
			(perimeter - Math.min(perimeter1, perimeter2)) / Math.abs(perimeter2 - perimeter1),
		0), 1);

		let numPoints, pointFraction;
		if (this.numPoints1 === this.numPoints2) {
			numPoints = this.numPoints1;
			pointFraction = 0;
		} else {
			// Greater perimeter needs to map to a larger number of points.
			if (this.numPoints2 < this.numPoints1) {
				t = 1 - t;
			}
			numPoints = this.numPoints1 * (1 - t) + this.numPoints2 * t;
			pointFraction = numPoints - Math.trunc(numPoints);
			numPoints = Math.trunc(numPoints);
		}

		this.lines = [];
		let increment = ((numPoints + pointFraction) / this.numPoints1 * this.increment) % numPoints;
		if (increment === 0) {
			return;
		}
		const mod1 = Math.ceil(this.numPoints1 / this.increment) * this.increment - this.numPoints1;
		const steps = Math.ceil(numPoints / increment);
		increment += (mod1 - (steps * increment - numPoints)) / steps;

		let offset = this.offset1 * (1 - t) + this.offset2 * t;
		offset += 1 / (numPoints + pointFraction);

		const epilson = 1 / 8192;
		increment = Math.round(increment / epilson) * epilson;
		if (increment === Math.trunc(increment)) {
			const divisor = Geometry.gcd(increment, numPoints);
			numPoints /= divisor;
			increment /= divisor;
		}

		const vertexFractions = Geometry.getPerimeterOffsets(
			morph.pointsX, morph.pointsY, this.closed
		);

		let colours, colourIndex, colour, prevColour, colourElapsed, colourLength;
		if (this.colourMorphs) {
			colours = this.colourMorphs.map(m => m.interpolate(morph, interpolation));
			colourIndex = 0;
			colour = colours[0];
			prevColour = colour;
			colourElapsed = 0;
			colourLength = this.colourLengths[0];
		}

		let numPieces = 1;
		if (numPoints % increment === 0) {
			numPieces = increment > 0.5 * numPoints ? numPoints - increment : increment;
		}

		for (let j = 0; j < numPieces; j++) {
			let start = j;
			const begin = start;
			for (let i = 0; i < numPoints / numPieces; i++) {
				const end = start + increment;
				if (end + offset > numPoints && !this.closed) {
					break;
				}

				const [vertex1, proportion1, x1, y1] = Geometry.getPerimeterOffset(
					((start + offset) % numPoints) / numPoints,
					morph.pointsX, morph.pointsY, vertexFractions, this.closed
				);
				const [vertex2, proportion2, x2, y2] = Geometry.getPerimeterOffset(
					((end + offset) % numPoints) / numPoints,
					morph.pointsX, morph.pointsY, vertexFractions, this.closed
				);

				if (colours) {
					const lineLength = Math.hypot(x2 - x1, y2 - y1);
					const lineColours = [colour];
					const colourOffsets = [0];
					let distance = colourLength - colourElapsed;
					let prevDistance = 0;
					while (distance < lineLength) {
						colourIndex = (colourIndex + 1) % colours.length;
						prevColour = colour;
						colour = colours[colourIndex];
						lineColours.push(colour);
						colourOffsets.push(distance / lineLength);
						colourLength = this.colourLengths[colourIndex];
						prevDistance = distance;
						distance += colourLength;
					}
					const percentage = (lineLength - prevDistance) / colourLength * 100;
					colour = 'color-mix(in srgb-linear, ' + prevColour + ', ' + colour + ' ' + percentage + '%)';
					lineColours.push(colour);
					colourOffsets.push(1);
					this.lines.push(new Line(x1, y1, x2, y2, lineColours, colourOffsets));
					colourElapsed = colourLength - (distance - lineLength);
				} else {
					this.lines.push(new Line(x1, y1, x2, y2));
				}

				start = end;
			}	// end for each length of string
		}	// end for each subshape

		if (this.startWidth !== undefined) {
			t = this.widthEase(interpolation);
			this.width = this.startWidth * (1 - t) + t * this.endWidth;
		}
	}

}

/**Represents calculation of a morph independently from the API or library used to draw the
 * shapes.
 */
class Morph {

	static DEFAULT_MAX_ROTATION = Math.PI / 2;

	constructor(
		polygon1, polygon2, fillMorph, strokeMorph, blendMode = 'source-over',
		startBlur = 0, endBlur = startBlur, maxRotation = Morph.DEFAULT_MAX_ROTATION
	) {
		polygon2.resize(polygon1.size);
		polygon2.rotate(polygon1, maxRotation);

		this.polygon1 = polygon1;
		this.polygon2 = polygon2;
		this.numFrames = 64;
		this.interpolationStep = 2 ** -6;

		this.translateX = 0;
		this.translateY = 0;
		this.rotation = 0;
		this.scaleX = 1;
		this.scaleY = 1;
		this.translateXEase = Ease.quadInOut;
		this.translateYEase = Ease.quadInOut;
		this.rotationEase = Ease.quadInOut;
		this.scaleXEase = Ease.linear;
		this.scaleYEase = Ease.linear;
		this.offsetEase = Ease.quadInOut;

		this.pointsX = undefined;
		this.pointsY = undefined;
		this.fillMorph = fillMorph;
		this.fillStyle = 'black';
		this.strokeMorph = strokeMorph;
		this.stroke = undefined;
		this.shadowMorph = undefined;
		// E.g. for producing string art.
		this.customMorph = undefined;

		this.blendMode = blendMode;
		this.startBlur = startBlur;
		this.endBlur = endBlur;
		this.blurEase = Ease.linear;
		this.blur = startBlur;
	}

	get numPoints() {
		return this.polygon2.numPoints;
	}

	setSpeed(speed) {
		// Technically the actual number of frames is one more than this.
		const numFrames = Math.ceil(moveLength(this.polygon1, this.polygon2) / speed);
		this.numFrames = numFrames;
		this.interpolationStep = 1 / numFrames;


		function noOp() {
			// Do nothing.
		}

		const numPoints = polygon2.numPoints;
		let hasPrecalculation = false;
		let fillState, strokeState;
		let fillPrecalc = noOp, strokePrecalc = noOp, fillEndPrecalc;
		let endFillPrecalc = noOp, endStrokePrecalc = noOp;
		if (this.fillMorph?.precalculate) {
			fillState = this.fillMorph.beginPrecalculation(numPoints);
			fillPrecalc = this.fillMorph.precalculate.bind(this.fillMorph);
			endFillPrecalc = this.fillMorph.endPrecalculation(this.fillMorph);
			hasPrecalculation = true;
		}
		if (this.strokeMorph?.colourMorph?.precalculate) {
			strokeState = this.strokeMorph.colourMorph.beginPrecalculation(numPoints)
			strokePrecalc = this.strokeMorph.colourMorph.precalculate.bind(this.strokeMorph.colourMorph);
			endStrokePrecalc = this.strokeMorph.colourMorph.endPrecalculation.bind(this.strokeMorph.colourMorph);
			hasPrecalculation = true;
		}
		if (hasPrecalculation) {
			const pointsX = new Array(numPoints);
			const pointsY = new Array(numPoints);
			for (let i = 0; i <= numFrames; i++) {
				const interpolation = i / numFrames;
				for (let j = 0; j < numPoints; j++) {
					pointsX[j] = polygon2.rotatedX[j] - (1 - interpolation) * polygon2.offsetsX[j];
					pointsY[j] = polygon2.rotatedY[j] - (1 - interpolation) * polygon2.offsetsY[j];
				}
				fillPrecalc(fillState, pointsX, pointsY, interpolation);
				strokePrecalc(strokeState, pointsX, pointsY, interpolation);
			}
			endFillPrecalc(fillState);
			endStrokePrecalc(strokeState);
		}
	}

	/**
	 * context {Object} Any object supporting createConicGradient, createLinearGradient and createRadialGradient.
	 */
	interpolate(interpolation, context, lightSourceX, lightSourceY) {
		const polygon1 = this.polygon1;
		const polygon2 = this.polygon2;
		let t;
		t = this.translateXEase(interpolation);
		this.translateX = polygon2.centreX * t + polygon1.centreX * (1 - t);
		t = this.translateYEase(interpolation);
		this.translateY = polygon2.centreY * t + polygon1.centreY * (1 - t);
		t = this.rotationEase(interpolation);
		this.rotation = -polygon2.rotation * t;
		t = this.scaleXEase(interpolation);
		this.scaleX = polygon2.scaleFactor * t + 1 - t;
		t = this.scaleYEase(interpolation);
		this.scaleY = polygon2.scaleFactor * t + 1 - t;

		const numPoints = polygon2.numPoints;
		const pointsX = new Array(numPoints);
		const pointsY = new Array(numPoints);
		t = this.offsetEase(interpolation);
		for (let i = 0; i < numPoints; i++) {
			pointsX[i] = polygon2.rotatedX[i] - (1 - t) * polygon2.offsetsX[i];
			pointsY[i] = polygon2.rotatedY[i] - (1 - t) * polygon2.offsetsY[i];
		}
		this.pointsX = pointsX;
		this.pointsY = pointsY;

		if (this.fillMorph) {
			this.fillStyle = this.fillMorph.interpolate(this, interpolation, context);
		}
		if (this.strokeMorph) {
			this.strokeMorph.interpolate(this, interpolation, context);
		}
		if (this.shadowMorph) {
			this.shadowMorph.interpolate(this, interpolation, lightSourceX, lightSourceY);
		}
		if (this.customMorph) {
			this.customMorph.interpolate(this, interpolation, context);
		}

		if (this.startBlur !== undefined) {
			t = this.blurEase(interpolation);
			this.blur = this.startBlur * (1 - t) + this.endBlur * t;
		}
	}

}

class Canvas2DMorph extends Morph {

	constructor(
		polygon1, polygon2, fillMorph, strokeMorph, blendMode = 'source-over',
		startBlur = 0, endBlur = startBlur, maxRotation = Morph.DEFAULT_MAX_ROTATION
	) {
		super(polygon1, polygon2, fillMorph, strokeMorph, blendMode, startBlur, endBlur, maxRotation);
	}

	transform(context) {
		context.translate(this.translateX, this.translateY);
		context.rotate(this.rotation);
		context.scale(this.scaleX, this.scaleY);

		context.fillStyle = this.fillStyle;
		const strokeMorph = this.strokeMorph;
		if (strokeMorph) {
			if (strokeMorph.colour !== undefined) {
				context.strokeStyle = strokeMorph.colour;
			}
			if (strokeMorph.width !== undefined) {
				context.lineWidth = strokeMorph.width;
			}
			if (strokeMorph.dash !== undefined) {
				context.setLineDash(strokeMorph.dash);
				context.lineDashOffset = strokeMorph.dashOffset;
			}
		}

		if (this.shadowMorph) {
			context.shadowOffsetX = this.shadowMorph.offsetX;
			context.shadowOffsetY = this.shadowMorph.offsetY;
			context.shadowBlur = this.shadowMorph.blur;
			context.shadowColor = this.shadowMorph.colour;
		}

		context.globalCompositeOperation = this.blendMode;
		if (this.blur !== undefined) {
			context.filter = 'blur(' + this.blur + 'px)';
		}
	}

	animate(
		context, mainCallback, beforeCallback = undefined, lightSourceX = 0, lightSourceY = 0
	) {
		const me = this;
		let startTime;

		function drawFrame(time) {
			let interpolation;
			if (startTime === undefined) {
				startTime = time;
				interpolation = 0;
			} else {
				// 60fps plus conversion from milliseconds.
				const frameNumber = (time - startTime) * 0.06;
				interpolation = Math.min(frameNumber * me.interpolationStep, 1);
			}
			me.interpolate(interpolation, context, lightSourceX, lightSourceY);

			const canvas = context.canvas;
			context.clearRect(0, 0, canvas.width, canvas.height);
			if (beforeCallback) {
				beforeCallback(context, me, interpolation);
			}
			context.save();
			me.transform(context);
			mainCallback(context, me, interpolation);
			context.restore();

			if (interpolation < 1) {
				requestAnimationFrame(drawFrame);
			}
		}
		requestAnimationFrame(drawFrame);
	}

	overlay(context, callback, lightSourceX = 0, lightSourceY = 0) {
		for (
			let interpolation = 0;
			interpolation <= 1;
			interpolation += this.interpolationStep
		) {
			this.interpolate(interpolation, context, lightSourceX, lightSourceY);
			context.save();
			this.transform(context);
			callback(context, this, interpolation);
			context.restore();
		}
	}

}

function drawSourceShapes(context, morph, interpolation) {
	const polygon1 = morph.polygon1;
	const polygon2 = morph.polygon2;

	Canvas.tracePolygonPath(context, polygon1.pointsX, polygon1.pointsY);
	context.fillStyle = 'rgba(255, 255, 192, 70%)';
	context.fill();

	Canvas.tracePolygonPath(context, polygon2.pointsX, polygon2.pointsY);
	context.fillStyle = 'rgba(224, 192, 192, 39%)';
	context.fill();

	if (morph.shadowMorph) {
		context.beginPath();
		context.arc(sunX, sunY, 35, 0, 2 * Math.PI);
		context.fillStyle = 'rgba(255, 192, 0, 30%)';
		context.fill();
	}
}

function drawShapeAndStrings(context, morph) {
	Canvas.drawInterpolatedShape(context, morph);
	Canvas.drawStrings(context, morph);
}

function drawFaded(context, morph, interpolation) {
	const alpha = Math.max(Math.round(255 / morph.numFrames), 1) / 255;
	context.globalAlpha = alpha;
	Canvas.drawInterpolatedShape(context, morph, interpolation);
	if (morph.customMorph) {
		Canvas.drawStrings(context, morph);
	}
}

const canvas = document.getElementById('canvas');
const context = canvas.getContext('2d');
canvas.width = window.innerWidth;
canvas.height = window.innerHeight;
context.lineJoin = 'round';

const Mode = Object.freeze({
	ANIMATE: 0,
	OVERLAY: 1,
});

const parameters = new URLSearchParams(document.location.search);
let speed = parseFloat(parameters.get('speed'));
let mode, fillMorph;

const fillStr = parameters.get('fill');
if (fillStr) {
	// Use %25 to write % inside a URL, %2B to write +, %23 to write #.
	fillMorph = new ColourMorph.Colour(fillStr);
}

switch (parameters.get('render')) {
case 'overlay':
	mode = Mode.OVERLAY;
	speed ||= 50;
	if (!fillMorph) {
		fillMorph = new ColourMorph.Colour(
			'rgb(calc(255-510*x), min(255*x, 255-255*x), calc(510*x-255))'
		);
	}
	break;
default:
	mode = Mode.ANIMATE;
	speed ||= 5;
	if (!fillMorph) {
		fillMorph = new ColourMorph.Constant('#00ff0063');
	}
}

const blendMode = parameters.get('blend') || 'source-over';

let startBlur = parseFloat(parameters.get('blur'));
let endBlur = parseFloat(parameters.get('blur2'));
if (!Number.isFinite(startBlur)) {
	startBlur = undefined;
}
if (!Number.isFinite(endBlur)) {
	endBlur = startBlur;
}

const numVertices = parseInt(parameters.get('vertices')) || 3 + Math.trunc(Math.random() * 6);
const numVertices2 = parseInt(parameters.get('vertices2')) || 3 + Math.trunc(Math.random() * 6);

let [polygon1, polygon2] = layoutPolygons(RandomShape.cyclicPolygonPoints, numVertices, numVertices2);
if (mode === Mode.OVERLAY && polygon2.scaleFactor > 1) {
	[polygon1, polygon2] = [polygon2, polygon1];
}

const fillStr2 = parameters.get('fill2');
if (fillStr2) {
	const toColour = new ColourMorph.Colour(fillStr2);
	const gradientType = parseInt(parameters.get('gradient'));
	/* 1 Edge to vertex
	 * 2 Along an edge
	 * 3 Vertex to vertex
	 * 4 Around a vertex
	 * 5 Around the centre
	 */
	switch (gradientType) {
	case 2:
		fillMorph = new ColourMorph.EdgeParallelGradient(0, [fillMorph, toColour]);
		break;
	case 3:
		const toVertex = Math.trunc(numVertices / 2);
		fillMorph = new ColourMorph.VertexVertexGradient(0, toVertex, [fillMorph, toColour]);
		break;
	case 4:
		fillMorph = new ColourMorph.VertexConicGradient(0, [fillMorph, toColour]);
		break;
	case 5:
		const spin = (parseFloat(parameters.get('spin')) || 0) * Math.PI / 180;
		fillMorph = new ColourMorph.CentreConicGradient(
			[fillMorph, toColour, fillMorph], [0, 0.5, 1], 0, spin
		);
		fillMorph.angleEase = Ease.expOut(3);
		break;
	default: // aka Case 1
		fillMorph = new ColourMorph.EdgePerpendicularGradient(0, [fillMorph, toColour]);
	}
}

const DEFAULT_LINE_WIDTH = 2;
const hasStrokeColour = parameters.has('stroke');
let startLineWidth = parseFloat(parameters.get('line_width'));
let endLineWidth = parseFloat(parameters.get('line_width2'));
const hasStartLineWidth = startLineWidth >= 0;
const hasEndLineWidth = endLineWidth >= 0;
if (!Number.isFinite(startLineWidth)) {
	startLineWidth = undefined;
}
if (!Number.isFinite(endLineWidth)) {
	endLineWidth = undefined;
}

let startDashPattern = [], endDashPattern = [];
let startDashOffset = 0, endDashOffset = 0;
if (parameters.has('dash')) {
	startDashPattern = parameters.get('dash').split(',').map(x => parseFloat(x));
	startDashOffset = parameters.get('dash_offset') || 0;
	endDashPattern = startDashPattern;
	endDashOffset = parameters.get('dash_offset2');
	if (endDashOffset === null) {
		endDashOffset = startDashOffset;
	}
}
if (parameters.has('dash2')) {
	endDashPattern = parameters.get('dash2').split(',').map(x => parseFloat(x));
}

let startStartStroke = parseFloat(parameters.get('start'));
if (!Number.isFinite(startStartStroke)) {
	startStartStroke = 0;
}
let endStartStroke = parseFloat(parameters.get('start2'));
if (!Number.isFinite(endStartStroke)) {
	endStartStroke = startStartStroke === 1 ? 0 : startStartStroke
}
let startEndStroke = parseFloat(parameters.get('end'));
if (!Number.isFinite(startEndStroke)) {
	startEndStroke = 1;
}
let endEndStroke = parseFloat(parameters.get('end2'))
if (!Number.isFinite(endEndStroke)) {
	endEndStroke = startEndStroke === 0 ? 1 : startEndStroke;
}

let strokeColourMorph, strokeColourMorph2, strokePaint;
let strokeMorph;

if (hasEndLineWidth && !hasStartLineWidth) {
	startLineWidth = DEFAULT_LINE_WIDTH;
}
if (hasStrokeColour) {
	strokeColourMorph = new ColourMorph.Colour(parameters.get('stroke'));
	strokeColourMorph2 = strokeColourMorph;
	strokePaint = strokeColourMorph;

	if (parameters.has('stroke2')) {
		strokeColourMorph2 = new ColourMorph.Colour(parameters.get('stroke2'));
		const gradientType = parseInt(parameters.get('stroke_gradient'));
		switch (gradientType) {
		case 1:
			strokePaint = new ColourMorph.EdgePerpendicularGradient(
				0, [strokeColourMorph, strokeColourMorph2]
			);
			break;
		case 2:
			strokePaint = new ColourMorph.EdgeParallelGradient(
				0, [strokeColourMorph, strokeColourMorph2]
			);
			break;
		default: // aka Case 5
			strokePaint = new ColourMorph.CentreConicGradient(
				[strokeColourMorph, strokeColourMorph2, strokeColourMorph],
				[0, 0.5, 1]
			);
		}
	}
}
if (
	!hasStartLineWidth && !hasEndLineWidth && !parameters.has('strings') &&
	(
		hasStrokeColour || startDashPattern.length > 0 ||
		startStartStroke > 0 || endStartStroke > 0 || startEndStroke < 1 || endEndStroke < 1
	)
) {
	startLineWidth = DEFAULT_LINE_WIDTH;
}

if (startLineWidth > 0 || endLineWidth > 0) {
	strokeMorph = new StrokeMorph(true, strokePaint);
	strokeMorph.setWidth(startLineWidth, endLineWidth);
	strokeMorph.setDash(startDashPattern, endDashPattern);
	strokeMorph.setDashOffset(startDashOffset, endDashOffset);
	strokeMorph.setStart(startStartStroke, endStartStroke);
	strokeMorph.setEnd(startEndStroke, endEndStroke);
}

let maxRotation = parseInt(parameters.get('max_rotation'));
if (Number.isFinite(maxRotation)) {
	maxRotation *= Math.PI / 180;
} else {
	maxRotation = Morph.DEFAULT_MAX_ROTATION;
}

let startShadowLength = parseFloat(parameters.get('shadow'));
if (!Number.isFinite(startShadowLength)) {
	startShadowLength = 0;
}

let endShadowLength = parseFloat(parameters.get('shadow2'));
if (!Number.isFinite(endShadowLength)) {
	endShadowLength = startShadowLength;
}

let shadowBlur = parseFloat(parameters.get('shadow_blur'));
if (!Number.isFinite(shadowBlur)) {
	shadowBlur = 0.4;
}

let shadowMorph = undefined;
if (startShadowLength !== 0 || endShadowLength !== 0) {
	const shadowColourStr = parameters.get('shadow_colour');
	const shadowColourMorph = shadowColourStr ?
		new ColourMorph.Colour(shadowColourStr) : ColourMorph.Constant.BLACK;
	shadowMorph = new ShadowMorph(
		startShadowLength, endShadowLength, shadowBlur, shadowColourMorph
	);
}

let sunX = parseFloat(parameters.get('sun_x'));
if (!Number.isFinite(sunX)) {
	sunX = 0.5;
}
sunX *= canvas.width;
let sunY = parseFloat(parameters.get('sun_y')) || 0;
sunY *= canvas.height;

const morph = new Canvas2DMorph(
	polygon1, polygon2, fillMorph, strokeMorph, blendMode, startBlur, endBlur,
	maxRotation
);
morph.shadowMorph = shadowMorph;

const startNumStringPoints = parseInt(parameters.get('strings'));
if (Number.isFinite(startNumStringPoints)) {
	let endNumStringPoints = parseInt(parameters.get('strings2'));
	if (!Number.isFinite(endNumStringPoints)) {
		endNumStringPoints = startNumStringPoints;
	}
	let increment = parseFloat(parameters.get('string_add'));
	if (!Number.isFinite(increment)) {
		increment = Math.round(startNumStringPoints / numVertices);
		increment += (1 - startNumStringPoints % increment) /
			Math.trunc(startNumStringPoints / increment);
	}
	const stringArt = new StringArtMorph(
		startNumStringPoints, 0, increment, true, endNumStringPoints
	);
	stringArt.setColours(
		[
			strokeColourMorph || ColourMorph.Constant.BLACK,
			strokeColourMorph2 || ColourMorph.Constant.BLACK
		],
		[200, 200]
	);
	morph.customMorph = stringArt;
}

morph.setSpeed(speed);

let translateEase = Ease.quadInOut;
morph.translateXEase = translateEase;
morph.translateYEase = translateEase;

switch (mode) {
case Mode.ANIMATE:
	// Draw as an animation.
	const drawFunction = morph.customMorph ? drawShapeAndStrings : Canvas.drawInterpolatedShape;
	morph.animate(context, drawFunction, drawSourceShapes, sunX, sunY);
	break;
case Mode.OVERLAY:
	// Draw with successive frames overlaid on top of each other.
	morph.overlay(context, drawFaded, sunX, sunY);
	break;
}

window.drawVertexMap = Canvas.drawVertexMap.bind(window, context, polygon1, polygon2);
