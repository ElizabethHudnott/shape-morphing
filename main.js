import {Geometry} from './src/math.js';
import {Shape, RandomShape} from './src/shapes.js';
import Ease from './src/easing.js';

const DEFAULT_MAX_ROTATION = Math.PI / 2;

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
	return Math.random() < 0.5 ? [shape1, shape2] : [shape2, shape1];
}

function polygonPath(context, pointsX, pointsY) {
	const numPoints = pointsX.length;
	context.beginPath();
	context.moveTo(pointsX[0], pointsY[0]);
	for (let i = 1; i < numPoints; i++) {
		const x = pointsX[i];
		const y = pointsY[i];
		context.lineTo(x, y);
	}
}

function polylinePath(context, pointsX, pointsY, closed = true, start = 0, end = 1) {
	const numPoints = pointsX.length;
	const segmentOffsets = new Array(numPoints + closed);
	segmentOffsets[0] = 0;
	let totalLength = 0;
	let prevX = pointsX[0];
	let prevY = pointsY[0];
	for (let i = 1; i < numPoints; i++) {
		const x = pointsX[i];
		const y = pointsY[i];
		totalLength += Math.hypot(x - prevX, y - prevY);
		segmentOffsets[i] = totalLength;
		prevX = x;
		prevY = y;
	}
	if (closed) {
		totalLength += Math.hypot(pointsX[0] - prevX, pointsY[0] - prevY);
		segmentOffsets[numPoints] = totalLength;
	}
	const startOffset = start * totalLength;
	let startVertex = 0;
	while (
		startVertex < numPoints + closed - 2 &&
		startOffset > segmentOffsets[startVertex + 1]
	) {
		startVertex++;
	}
	const endOffset = end * totalLength;
	let endVertex = startVertex;
	while (
		endVertex < numPoints + closed - 2 &&
		endOffset > segmentOffsets[endVertex + 1]
	) {
		endVertex++;
	}

	const startX1 = pointsX[startVertex];
	const startY1 = pointsY[startVertex];
	let nextVertex = startVertex + 1;
	let add = 0;
	if (nextVertex === numPoints) {
		nextVertex = 0;
		add = totalLength;
	}
	const startX2 = pointsX[nextVertex];
	const startY2 = pointsY[nextVertex];
	const startOffset1 = segmentOffsets[startVertex];
	const startOffset2 = segmentOffsets[nextVertex] + add;
	const startProportion = (startOffset - startOffset1) / (startOffset2 - startOffset1);
	const startX = startX1 * (1 - startProportion) + startX2 * startProportion;
	const startY = startY1 * (1 - startProportion) + startY2 * startProportion;

	const endX1 = pointsX[endVertex];
	const endY1 = pointsY[endVertex];
	nextVertex = endVertex + 1;
	add = 0;
	if (nextVertex === numPoints) {
		nextVertex = 0;
		add = totalLength;
	}
	const endX2 = pointsX[nextVertex];
	const endY2 = pointsY[nextVertex];
	const endOffset1 = segmentOffsets[endVertex];
	const endOffset2 = segmentOffsets[nextVertex] + add;
	const endProportion = (endOffset - endOffset1) / (endOffset2 - endOffset1);
	const endX = endX1 * (1 - endProportion) + endX2 * endProportion;
	const endY = endY1 * (1 - endProportion) + endY2 * endProportion;

	context.beginPath();
	context.moveTo(startX, startY);
	for (let i = startVertex + 1; i <= endVertex; i++) {
		context.lineTo(pointsX[i], pointsY[i]);
	}
	context.lineTo(endX, endY);
	if (closed && start === 0 && end === 1) {
		context.closePath();
	}
}

function moveLength(polygon1, polygon2) {
	const numPoints = polygon2.numPoints;
	const inverseScale = polygon1.size / polygon2.size;

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

class ConstantColour {

	static BLACK = new ConstantColour('black');

	constructor(colour) {
		this.colour = colour;
	}

	interpolate(morph, interpolation) {
		return this.colour;
	}
}

class ColourMorph {

	constructor(str, easings = (new Array(4)).fill(Ease.linear)) {
		this.str = str.replace(/\s/g, '').replace(/[+\-]/g, ' $& ');
		this.easings = easings;
	}

	interpolate(morph, interpolation) {
		const values = this.easings.map(f => f(interpolation));
		let str = this.str.replace(/\bw\b/g, values[0]);
		str = str.replace(/\bx\b/g, values[1]);
		str = str.replace(/\by\b/g, values[2]);
		str = str.replace(/\bz\b/g, values[3]);
		return str;
	}

}

class VertexGradientMorph {

	constructor(vertexNum1, vertexNum2, colourMorphs, offsets = [0, 1]) {
		this.vertexNum1 = vertexNum1;
		this.vertexNum2 = vertexNum2;
		this.colourMorphs = colourMorphs;
		this.offsets = offsets;
	}

	interpolate(morph, interpolation, context) {
		const x1 = morph.pointsX[this.vertexNum1];
		const y1 = morph.pointsY[this.vertexNum1];
		const x2 = morph.pointsX[this.vertexNum2];
		const y2 = morph.pointsY[this.vertexNum2];
		const gradient = context.createLinearGradient(x1, y1, x2, y2);
		for (let i = 0; i < this.offsets.length; i++) {
			const colour = this.colourMorphs[i].interpolate(morph, interpolation);
			gradient.addColorStop(this.offsets[i], colour);
		}
		return gradient;
	}

}

class EdgeParallelGradientMorph {

	constructor(vertexNum, colourMorphs, offsets = [0, 1]) {
		this.vertexNum = vertexNum;
		this.colourMorphs = colourMorphs;
		this.offsets = offsets;
	}

	interpolate(morph, interpolation, context) {
		const pointsX = morph.pointsX;
		const pointsY = morph.pointsY;
		const numPoints = pointsX.length;
		const vertex1 =  this.vertexNum;
		const vertex2 = (vertex1 + 1) % numPoints;
		const edgeX1 = pointsX[vertex1];
		const edgeY1 = pointsY[vertex1];
		const edgeX2 = pointsX[vertex2];
		const edgeY2 = pointsY[vertex2];
		const [edgeDeltaX, a1, b1, c1, a2, b2] = Geometry.projectionOntoLine(edgeX1, edgeY1, edgeX2, edgeY2);

		// Min and max x values and associated y values or min and max y and associated
		// constant x value.
		let minX, maxX, minY, maxY;
		let invertDirection = false;
		if (edgeX1 < edgeX2 || (edgeDeltaX === 0 && edgeY1 < edgeY2)) {
			[minX, minY, maxX, maxY] = [edgeX1, edgeY1, edgeX2, edgeY2];
		} else {
			[minX, minY, maxX, maxY] = [edgeX2, edgeY2, edgeX1, edgeY1];
			invertDirection = true;
		}

		for (let i = 1; i <= numPoints - 2; i++) {
			const targetVertex = (vertex2 + i) % numPoints;
			const targetX = pointsX[targetVertex];
			const targetY = pointsY[targetVertex];
			const [intersectX, intersectY] = Geometry.projectOntoLine(targetX, targetY, a1, b1, c1, a2, b2);

			if (edgeDeltaX !== 0) {
				if (intersectX < minX) {
					[minX, minY] = [intersectX, intersectY];
				} else if (intersectX > maxX) {
					[maxX, maxY] = [intersectX, intersectY];
				}
			} else {
				if (intersectY < minY) {
					minY = intersectY;
				} else if (intersectY > maxY) {
					maxY = intersectY;
				}
			}

		}

		let gradient;
		if (!invertDirection) {
			gradient = context.createLinearGradient(minX, minY, maxX, maxY);
		} else {
			gradient = context.createLinearGradient(maxX, maxY, minX, minY);
		}
		for (let i = 0; i < this.offsets.length; i++) {
			const colour = this.colourMorphs[i].interpolate(morph, interpolation);
			gradient.addColorStop(this.offsets[i], colour);
		}
		return gradient;
	}

}

class EdgePerpendicularGradientMorph {

	constructor(vertexNum, colourMorphs, offsets = [0, 1]) {
		this.vertexNum = vertexNum;
		this.targetVertex = undefined;
		this.colourMorphs = colourMorphs;
		this.offsets = offsets;
	}

	beginPrecalculation(numPoints) {
		const counts = new Array(numPoints);
		counts.fill(0);
		return counts;
	}

	precalculate(counts, pointsX, pointsY, interpolation) {
		const numPoints = pointsX.length;
		const vertex1 =  this.vertexNum;
		const vertex2 = (vertex1 + 1) % numPoints;
		const edgeX1 = pointsX[vertex1];
		const edgeY1 = pointsY[vertex1];
		const edgeX2 = pointsX[vertex2];
		const edgeY2 = pointsY[vertex2];
		const [edgeDeltaX, a1, b1, c1, a2, b2] = Geometry.projectionOntoLine(edgeX1, edgeY1, edgeX2, edgeY2);

		let maxDistanceSq = 0;
		let targetVertex = 0;
		for (let i = 0; i < numPoints; i++) {
			const targetX = pointsX[i];
			const targetY = pointsY[i];
			const [intersectX, intersectY] = Geometry.projectOntoLine(targetX, targetY, a1, b1, c1, a2, b2);
			const [measureX, measureY] = Geometry.constrainToLineSegment(
				intersectX, intersectY, edgeX1, edgeY1, edgeX2, edgeY2
			);

			const gradientDeltaX = targetX - measureX;
			const gradientDeltaY = targetY - measureY;
			const distanceSq = gradientDeltaX * gradientDeltaX + gradientDeltaY * gradientDeltaY;
			if (distanceSq > maxDistanceSq) {
				targetVertex = i;
				maxDistanceSq = distanceSq;
			}
		}
		counts[targetVertex]++;
	}

	endPrecalculation(counts) {
		const numPoints = counts.length;
		let targetVertex = 0;
		let maxCount = 0;
		for (let i = 0; i < numPoints; i++) {
			const count = counts[i];
			if (count > maxCount) {
				targetVertex = i;
				maxCount = count;
			}
		}
		this.targetVertex = targetVertex;
	}

	interpolate(morph, interpolation, context) {
		const pointsX = morph.pointsX;
		const pointsY = morph.pointsY;
		const numPoints = pointsX.length;
		const vertex1 =  this.vertexNum;
		const vertex2 = (vertex1 + 1) % numPoints;
		const edgeX1 = pointsX[vertex1];
		const edgeY1 = pointsY[vertex1];
		const edgeX2 = pointsX[vertex2];
		const edgeY2 = pointsY[vertex2];
		const [edgeDeltaX, a1, b1, c1, a2, b2] = Geometry.projectionOntoLine(edgeX1, edgeY1, edgeX2, edgeY2);

		const x2 = pointsX[this.targetVertex];
		const y2 = pointsY[this.targetVertex];
		const c2 = -a2 * x2 - b2 * y2;
		const denominator = a1 * b2 - a2 * b1;
		const x1 = (b1 * c2 - b2 * c1) / denominator;
		const y1 = (c1 * a2 - c2 * a1) / denominator;

		const gradient = context.createLinearGradient(x1, y1, x2, y2);
		for (let i = 0; i < this.offsets.length; i++) {
			const colour = this.colourMorphs[i].interpolate(morph, interpolation);
			gradient.addColorStop(this.offsets[i], colour);
		}
		return gradient;
	}

}

class VertexConicGradientMorph {

	constructor(vertexNum, colourMorphs, offsets = [0, 1]) {
		this.vertexNum = vertexNum;
		this.colourMorphs = colourMorphs;
		this.offsets = offsets;
	}

	interpolate(morph, interpolation, context) {
		const pointsX = morph.pointsX;
		const pointsY = morph.pointsY;
		const numPoints = pointsX.length
		const vertex1 = (this.vertexNum + numPoints - 1) % numPoints;
		const vertex2 =  this.vertexNum;
		const vertex3 = (vertex2 + 1) % numPoints;
		const x1 = pointsX[vertex1];
		const y1 = pointsY[vertex1];
		const x2 = pointsX[vertex2];
		const y2 = pointsY[vertex2];
		const x3 = pointsX[vertex3];
		const y3 = pointsY[vertex3];
		const startAngle = Math.atan2(y1 - y2, x1 - x2);
		let endAngle = Math.atan2(y3 - y2, x3 - x2);
		if (endAngle < startAngle) {
			endAngle += 2 * Math.PI;
		}
		const turnAmount = (endAngle - startAngle) / (2 * Math.PI);
		const gradient = context.createConicGradient(startAngle, x2, y2);
		for (let i = 0; i < this.offsets.length; i++) {
			const colour = this.colourMorphs[i].interpolate(morph, interpolation);
			gradient.addColorStop(turnAmount * this.offsets[i], colour);
		}
		return gradient;
	}

}

class CentreConicGradientMorph {

	constructor(colourMorphs, offsets = [0, 1], startAngle = 0, endAngle = startAngle) {
		this.startAngle = startAngle;
		this.endAngle = endAngle;
		this.angleEase = Ease.linear;
		this.colourMorphs = colourMorphs;
		this.offsets = offsets;
	}

	interpolate(morph, interpolation, context) {
		const t = this.angleEase(interpolation);
		const angle = this.endAngle * t + this.startAngle * (1 - t);
		const gradient = context.createConicGradient(angle, 0, 0);
		for (let i = 0; i < this.offsets.length; i++) {
			const colour = this.colourMorphs[i].interpolate(morph, interpolation);
			gradient.addColorStop(this.offsets[i], colour);
		}
		return gradient;
	}

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

	constructor(startLength, endLength, blur = 0.25, colourMorph = ConstantColour.BLACK) {
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

class Morph {

	constructor(
		polygon1, polygon2, fillMorph, strokeMorph, shadowMorph, blendMode = 'source-over',
		startBlur = 0, endBlur = startBlur, maxRotation = DEFAULT_MAX_ROTATION
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
		this.shadowMorph = shadowMorph;
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

		if (this.fillMorph?.precalculate) {
			const numPoints = polygon2.numPoints;
			const pointsX = new Array(numPoints);
			const pointsY = new Array(numPoints);
			const state = this.fillMorph.beginPrecalculation(numPoints);
			for (let i = 0; i <= numFrames; i++) {
				const interpolation = i / numFrames;
				for (let j = 0; j < numPoints; j++) {
					pointsX[j] = polygon2.rotatedX[j] - (1 - interpolation) * polygon2.offsetsX[j];
					pointsY[j] = polygon2.rotatedY[j] - (1 - interpolation) * polygon2.offsetsY[j];
				}
				this.fillMorph.precalculate(state, pointsX, pointsY, interpolation);
			}
			this.fillMorph.endPrecalculation(state);
		}
	}

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
		this.scaleX = polygon2.size / polygon1.size * t + 1 - t;
		t = this.scaleYEase(interpolation);
		this.scaleY = polygon2.size / polygon1.size * t + 1 - t;

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

		if (this.startBlur !== undefined) {
			t = this.blurEase(interpolation);
			this.blur = this.startBlur * (1 - t) + this.endBlur * t;
		}
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

	polygonPath(context, polygon1.pointsX, polygon1.pointsY);
	context.fillStyle = 'rgba(255, 255, 192, 70%)';
	context.fill();

	polygonPath(context, polygon2.pointsX, polygon2.pointsY);
	context.fillStyle = 'rgba(224, 192, 192, 39%)';
	context.fill();

	if (morph.shadowMorph) {
		context.beginPath();
		context.arc(sunX, sunY, 35, 0, 2 * Math.PI);
		context.fillStyle = 'rgba(255, 192, 0, 30%)';
		context.fill();
	}
}

function drawInterpolatedShape(context, morph, interpolation) {
	polygonPath(context, morph.pointsX, morph.pointsY);
	context.fill();

	const strokeMorph = morph.strokeMorph;
	if (strokeMorph && strokeMorph.width > 0) {
		polylinePath(
			context, morph.pointsX, morph.pointsY, strokeMorph.closed, strokeMorph.start,
			strokeMorph.end
		);
		context.globalAlpha = 1;
		context.shadowColor = 'transparent';
		context.stroke();
	}
}

function drawFaded(context, morph, interpolation) {
	const alpha = Math.max(Math.round(1 * 255 / morph.numFrames), 1) / 255;
	context.globalAlpha = alpha;
	drawInterpolatedShape(context, morph, interpolation);
}

function drawVertexMap() {
	const numPoints = polygon1.numPoints;
	context.clearRect(0, 0, canvas.width, canvas.height);
	context.save();
	context.translate(0.5 * canvas.width,  0.5 * canvas.height);

	context.beginPath();
	context.moveTo(polygon1.deltaX[0], polygon1.deltaY[0]);
	for (let i = 1; i < numPoints; i++) {
		context.lineTo(polygon1.deltaX[i], polygon1.deltaY[i]);
	}
	context.closePath();
	context.fillStyle = 'rgba(255, 0, 0, 70%)';
	context.fill();

	context.beginPath();
	context.moveTo(polygon2.rotatedX[0], polygon2.rotatedY[0]);
	for (let i = 1; i < numPoints; i++) {
		context.lineTo(polygon2.rotatedX[i], polygon2.rotatedY[i]);
	}
	context.closePath();
	context.fillStyle = 'rgba(0, 255, 0, 70%)';
	context.fill();

	for (let i = 0; i < numPoints; i++) {
		const x2 = polygon2.rotatedX[i];
		const y2 = polygon2.rotatedY[i]
		const x1 = x2 - polygon2.offsetsX[i];
		const y1 = y2 - polygon2.offsetsY[i];
		context.beginPath();
		context.moveTo(x1, y1);
		context.lineTo(x2, y2);
		const gradient = context.createLinearGradient(x1, y1, x2, y2);
		gradient.addColorStop(0, 'black');
		gradient.addColorStop(1, 'blue');
		context.strokeStyle = gradient;
		context.stroke();
	}
	context.restore();
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
	fillMorph = new ColourMorph(fillStr);
}

switch (parameters.get('render')) {
case 'overlay':
	mode = Mode.OVERLAY;
	speed ||= 50;
	if (!fillMorph) {
		fillMorph = new ColourMorph('rgb(calc(255-510*x), min(255*x, 255-255*x), calc(510*x-255))');
	}
	break;
default:
	mode = Mode.ANIMATE;
	speed ||= 5;
	if (!fillMorph) {
		fillMorph = new ConstantColour('#00ff0063');
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
if (mode === Mode.OVERLAY && polygon2.size > polygon1.size) {
	[polygon1, polygon2] = [polygon2, polygon1];
}

const fillStr2 = parameters.get('fill2');
if (fillStr2) {
	const toColour = new ColourMorph(fillStr2);
	const gradientType = parseInt(parameters.get('gradient'));
	/* 1 Edge to vertex
	 * 2 Along an edge
	 * 3 Vertex to vertex
	 * 4 Around a vertex
	 * 5 Around the centre
	 */
	switch (gradientType) {
	case 2:
		fillMorph = new EdgeParallelGradientMorph(0, [fillMorph, toColour]);
		break;
	case 3:
		const toVertex = Math.trunc(numVertices / 2);
		fillMorph = new VertexGradientMorph(0, toVertex, [fillMorph, toColour]);
		break;
	case 4:
		fillMorph = new VertexConicGradientMorph(0, [fillMorph, toColour]);
		break;
	case 5:
		const spin = (parseFloat(parameters.get('spin')) || 0) * Math.PI / 180;
		fillMorph = new CentreConicGradientMorph(
			[fillMorph, toColour, fillMorph], [0, 0.5, 1], 0, spin
		);
		fillMorph.angleEase = Ease.expOut(3);
		break;
	default:
		fillMorph = new EdgePerpendicularGradientMorph(0, [fillMorph, toColour]);
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

let strokeColourMorph;
let strokeMorph;

if (hasEndLineWidth && !hasStartLineWidth) {
	startLineWidth = DEFAULT_LINE_WIDTH;
}
if (hasStrokeColour) {
	strokeColourMorph = new ColourMorph(parameters.get('stroke'));
	if (parameters.has('stroke2')) {
		toColour = new ColourMorph(parameters.get('stroke2'));
		const gradientType = parseInt(parameters.get('stroke_gradient'));
		switch (gradientType) {
		case 2:
			strokeColourMorph = new EdgeParallelGradientMorph(0, [strokeColourMorph, toColour]);
			break;
		default:
			strokeColourMorph = new CentreConicGradientMorph(
				[strokeColourMorph, toColour, strokeColourMorph],
				[0, 0.5, 1]
			);
		}
	}
}
if (
	(!hasStartLineWidth && !hasEndLineWidth) &&
	(
		hasStrokeColour || startDashPattern.length > 0 ||
		startStartStroke > 0 || endStartStroke > 0 || startEndStroke < 1 || endEndStroke < 1
	)
) {
	startLineWidth = DEFAULT_LINE_WIDTH;
}

if (startLineWidth > 0 || endLineWidth > 0) {
	strokeMorph = new StrokeMorph(true, strokeColourMorph);
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
	maxRotation = DEFAULT_MAX_ROTATION;
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
		new ColourMorph(shadowColourStr) : ConstantColour.BLACK;
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

const morph = new Morph(
	polygon1, polygon2, fillMorph, strokeMorph, shadowMorph, blendMode, startBlur, endBlur,
	maxRotation
);
morph.setSpeed(speed);

let translateEase = Ease.quadInOut;
morph.translateXEase = translateEase;
morph.translateYEase = translateEase;

switch (mode) {
case Mode.ANIMATE:
	// Draw as an animation.
	morph.animate(context, drawInterpolatedShape, drawSourceShapes, sunX, sunY);
	break;
case Mode.OVERLAY:
	// Draw with successive frames overlaid on top of each other.
	morph.overlay(context, drawFaded, sunX, sunY);
	break;
}

window.drawVertexMap = drawVertexMap;
