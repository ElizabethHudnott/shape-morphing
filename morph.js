const DEFAULT_MAX_ROTATION = Math.PI / 3;

const numericAscending = (a, b) => a - b;
const numericDescending = (a, b) => b - a;

function sortVertical(pointsX, pointsY, startIndex, endIndex) {
	let i = startIndex;
	let prevEndY;
	while (i < endIndex) {
		const startX = pointsX[i];
		const startY = pointsY[i];
		const subStartIndex = i;
		let endY = startY;
		let x;
		i++;
		do {
			const x = pointsX[i];
			const y = pointsY[i];
			if (x !== startX) {
				break;
			}
			if (y === endY) {
				pointsX.splice(i, 1);
				pointsY.splice(i, 1);
				endIndex--;
			} else {
				endY = y;
				i++;
			}
		} while (i <= endIndex);
		const subPointsY = pointsY.slice(subStartIndex, i);
		if (endY > prevEndY) {
			subPointsY.sort(numericAscending);
		} else {
			subPointsY.sort(numericDescending);
		}
		pointsY.splice(subStartIndex, subPointsY.length, ...subPointsY);
		prevEndY = endY;
	}
}

function orderPolygonPoints(pointsX, pointsY) {
	let numPoints = pointsX.length;
	let left = pointsX[0];
	let right = left;
	let leftIndex = 0;
	let rightIndex = 0;
	for (let i = 1; i < numPoints; i++) {
		const x = pointsX[i];
		if (x < left) {
			left = x;
			leftIndex = i;
		} else if (x > right) {
			right = x;
			rightIndex = i;
		}
	}
	const leftY = pointsY[leftIndex];
	const rightY = pointsY[rightIndex];
	const gradient = (rightY - leftY) / (right - left);
	const intercept = leftY - gradient * left;

	const indexesAbove = [];
	const indexesBelow = [];
	for (let i = 0; i < numPoints; i++) {
		const x = pointsX[i];
		const y = pointsY[i];
		const lineY = gradient * x + intercept;
		if (y >= lineY) {
			indexesAbove.push(i);
		} else {
			indexesBelow.push(i);
		}
	}

	indexesAbove.sort( (i, j) => pointsX[i] - pointsX[j] );
	indexesBelow.sort( (i, j) => pointsX[j] - pointsX[i] );
	let numAbove = indexesAbove.length;
	const numBelow = indexesBelow.length;
	const newX = new Array(numPoints);
	const newY = new Array(numPoints);

	for (let i = 0; i < numAbove; i++) {
		const sourceIndex = indexesAbove[i];
		newX[i] = pointsX[sourceIndex];
		newY[i] = pointsY[sourceIndex];
	}
	sortVertical(newX, newY, 0, numAbove - 1);
	numPoints = newX.length;
	numAbove = numPoints - numBelow;

	for (let i = 0; i < numBelow; i++) {
		const sourceIndex = indexesBelow[i];
		const destIndex = i + numAbove;
		newX[destIndex] = pointsX[sourceIndex];
		newY[destIndex] = pointsY[sourceIndex];
	}
	sortVertical(newX, newY, numAbove, numPoints - 1);

	return [newX, newY];
}

function twiceTriangleArea(x1, y1, x2, y2, x3, y3) {
	return (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
}

function centroid(pointsX, pointsY) {
	const numPoints = pointsX.length;
	const x0 = pointsX[0];
	const y0 = pointsY[0];
	let doubleArea = 0, centreX = 0, centreY = 0;
	for (let i = 1; i < numPoints; i++) {
		const index2 = (i + 1) % numPoints;
		const x1 = pointsX[i];
		const y1 = pointsY[i];
		const x2 = pointsX[index2];
		const y2 = pointsY[index2];
		const doubleTriArea = twiceTriangleArea(x0, y0, x1, y1, x2, y2);
		centreX += doubleTriArea * (x0 + x1 + x2);
		centreY += doubleTriArea * (y0 + y1 + y2);
		doubleArea += doubleTriArea;
	}
	centreX /= 3 * doubleArea;
	centreY /= 3 * doubleArea;
	return [centreX, centreY];
}

class Shape {

	/**Constructs a new shape from two arrays of x and y coordinates and performs the following
	 * operations.
	 * 1) Finds the shape's centre. (this.centreX and this.centreY)
	 * 2) Finds the positions of the shape's vertices relative to its centre. (deltaX, deltaY)
	 * 3) Finds the polar representations of those displacements. (angles, radii)
	 * 4) Sorts the points by ascending order of their angle when expressed in polar form.
	 * 5) Computes the average radius. (size)
	 */
	constructor(pointsX, pointsY) {
		const numPoints = pointsX.length;
		[pointsX, pointsY] = orderPolygonPoints(pointsX, pointsY);
		const [centreX, centreY] = centroid(pointsX, pointsY);

		const deltaX = new Array(numPoints);
		const deltaY = new Array(numPoints);
		const angles = new Array(numPoints);
		const radii = new Array(numPoints);

		for (let i = 0; i < numPoints; i++) {
			const dx = pointsX[i] - centreX;
			const dy = pointsY[i] - centreY;
			deltaX[i] = dx;
			deltaY[i] = dy;
			radii[i] = Math.hypot(dx, dy);
			angles[i] = Math.atan2(deltaY[i], deltaX[i]);
		}

		let area = 0;
		for (let i = 0; i < numPoints; i++) {
			const angle1 = angles[i];
			const radius1 = radii[i];
			let angle2, radius2;
			if (i === numPoints - 1) {
				angle2 = angles[0] + 2 * Math.PI;
				radius2 = radii[0];
			} else {
				angle2 = angles[i + 1];
				radius2 = radii[i + 1];
			}
			// Area of a sector of an Archimedean spiral.
			// r = a + b * t
			// area = 0.5 * Integral r^2 dt [angle1, angle2]
			// r^2 = (a + bt) ^ 2 = a^2 + b^2 * t^2 + 2abt
			// Integral r^2 dt = a^2 * t + 1/3 * b^2 * t^3 + abt^2
			const angleDiff = angle2 - angle1;
			const b = (radius2 - radius1) / angleDiff;
			const a = radius1 - b * angle1;
			const aSquare = a * a;
			const bSquare = b * b;
			const angle1Square = angle1 * angle1;
			const angle2Square = angle2 * angle2;
			const angleSquareDiff = angle2Square - angle1Square;
			const angle1Cube = angle1Square * angle1;
			const angle2Cube = angle2Square * angle2;
			const angleCubeDiff = angle2Cube - angle1Cube;
			const sectorArea = aSquare * angleDiff + bSquare * angleCubeDiff / 3 + a * b * angleSquareDiff;
			area += Math.abs(sectorArea);
		}
		area *= 0.5;

		this.numPoints = numPoints;
		this.pointsX = pointsX;
		this.pointsY = pointsY;
		this.centreX = centreX;
		this.centreY = centreY;
		this.deltaX = deltaX;
		this.deltaY = deltaY;
		this.angles = angles;
		this.radii = radii;
		this.size = Math.sqrt(area / Math.PI);
		this.resetTransform();
	}

	resetTransform() {
		const deltaX = this.deltaX;
		const deltaY = this.deltaY;
		this.resizedX = deltaX;
		this.resizedY = deltaY;
		this.rotatedX = deltaX;
		this.rotatedY = deltaY;
		this.rotation = 0;
		// Align Vertex 0 of this shape with this numbered vertex in the target shape.
		this.vertexRotations = 0;
		// The offsets that remain after translation, resizing and rotating.
		this.offsetsX = undefined;
		this.offsetsY = undefined;
	}

	/**Computes a new set of points by scaling the shape's original points to match a different
	 * average radius and stores the rectangular coordinates of the new points in this.resizedX
	 * and this.resizedY
	 */
	resize(size) {
		const scale = size / this.size;
		const numPoints = this.numPoints;
		for (let i = 0; i < numPoints; i++) {
			const radius = this.radii[i] * scale;
			const angle = this.angles[i];
			this.resizedX[i] = radius * Math.cos(angle);
			this.resizedY[i] = radius * Math.sin(angle);
		}
	}

	/**Computes a new set of points by rotating the scaled shape to align with the target shape
	 * as best as possible and stores the rectangular coordinates of the new points in
	 * this.rotatedX and this.rotatedY, along with the residual displacement vectors in
	 * this.offsetsX and this.offsetsY.
	 */
	rotate(shape2, maxRotation) {
		const numPoints = this.numPoints;
		const numPoints2 = shape2.numPoints;
		let minError = Infinity;
		let minErrorVertex = 0;
		let selectedRotatedX = this.resizedX;
		let selectedRotatedY = this.resizedY;
		let selectedRotation = 0;
		for (let i = 0; i < numPoints2; i++) {
			let rotation = shape2.angles[i] - this.angles[0];
			if (rotation < -Math.PI) {
				rotation += 2 * Math.PI;
			} else if (rotation > Math.PI) {
				rotation -= 2 * Math.PI;
			}
			if (Math.abs(rotation) > maxRotation) {
				// Rotations bigger than a certain amount can be a bit dizzifying.
				continue;
			}

			const cos = Math.cos(rotation);
			const sin = Math.sin(rotation);
			const rotatedX = new Array(numPoints);
			const rotatedY = new Array(numPoints);
			const distancesSquared = new Array(numPoints);
			for (let j = 0; j < numPoints; j++) {
				const x = this.resizedX[j];
				const y = this.resizedY[j];
				const transformedX = x * cos - y * sin;
				const transformedY = x * sin + y * cos;
				const vertexIndex2 = (i + j) % numPoints2;
				const targetX = shape2.resizedX[vertexIndex2];
				const targetY = shape2.resizedY[vertexIndex2];
				const distanceX = targetX - transformedX;
				const distanceY = targetY - transformedY;
				distancesSquared[j] = distanceX * distanceX + distanceY * distanceY;
				rotatedX[j] = transformedX;
				rotatedY[j] = transformedY;
			}
			distancesSquared.sort(numericAscending);
			const midIndex = numPoints >> 1;
			let medianDistance;
			if (numPoints & 1) {
				medianDistance = distancesSquared[midIndex];
			} else {
				medianDistance = 0.5 * (distancesSquared[midIndex] + distancesSquared[midIndex - 1]);
			}
			if (medianDistance < minError) {
				minError = medianDistance;
				minErrorVertex = i;
				const x0 = this.resizedX[0];
				const y0 = this.resizedY[0];
				selectedRotatedX = rotatedX;
				selectedRotatedY = rotatedY;
				selectedRotatedX[0] = x0 * cos - y0 * sin;
				selectedRotatedY[0] = x0 * sin + y0 * cos;
				selectedRotation = rotation;
			}
		}
		this.rotatedX = selectedRotatedX;
		this.rotatedY = selectedRotatedY;
		this.rotation = selectedRotation;
		this.vertexRotations = minErrorVertex;
		const offsetsX = new Array(numPoints);
		const offsetsY = new Array(numPoints);
		for (let i = 0; i < numPoints; i++) {
			const vertexIndex2 = (i + minErrorVertex) % numPoints2;
			offsetsX[i] = selectedRotatedX[i] - shape2.resizedX[vertexIndex2];
			offsetsY[i] = selectedRotatedY[i] - shape2.resizedY[vertexIndex2];
		}
		this.offsetsX = offsetsX;
		this.offsetsY = offsetsY;
	}

}

function randomPolygon(numPoints) {
	const pointsX = new Array(numPoints);
	const pointsY = new Array(numPoints);
	for (let i = 0; i < numPoints; i++) {
		pointsX[i] = Math.trunc(Math.random() * canvas.width);
		pointsY[i] = Math.trunc(Math.random() * canvas.height);
	}
	return new Shape(pointsX, pointsY);
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
	constructor(colour) {
		this.colour = colour;
	}

	interpolate(morph, interpolation) {
		return this.colour;
	}
}

class ColourMorph {

	constructor(str) {
		this.str = str.replace(/\s/g, '').replace(/[+\-]/g, ' $& ');
	}

	interpolate(morph, interpolation) {
		return this.str.replace(/\bt\b/g, interpolation);
	}

}

class VertexGradientMorph {

	constructor(colourMorph1, colourMorph2, vertexNum1, vertexNum2) {
		this.colourMorph1 = colourMorph1;
		this.colourMorph2 = colourMorph2;
		this.vertexNum1 = vertexNum1;
		this.vertexNum2 = vertexNum2;
	}

	interpolate(morph, interpolation, context) {
		const x1 = morph.pointsX[this.vertexNum1];
		const y1 = morph.pointsY[this.vertexNum1];
		const x2 = morph.pointsX[this.vertexNum2];
		const y2 = morph.pointsY[this.vertexNum2];
		const colour1 = this.colourMorph1.interpolate(morph, interpolation);
		const colour2 = this.colourMorph2.interpolate(morph, interpolation);
		const gradient = context.createLinearGradient(x1, y2, x2, y2);
		gradient.addColorStop(0, colour1);
		gradient.addColorStop(1, colour2);
		return gradient;
	}

}

class EdgeGradientMorph {
	constructor(colourMorph1, colourMorph2, vertexNum) {
		this.colourMorph1 = colourMorph1;
		this.colourMorph2 = colourMorph2;
		this.vertexNum = vertexNum;
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
		const edgeDeltaX = edgeX2 - edgeX1;
		const edgeDeltaY = edgeY2 - edgeY1;
		// Line 1: Our edge
		// Line 2: Perpendicular to our edge running through another vertex
		// a1 * x + b1 * y + c1 = 0
		// a2 * x + b2 * y + c2 = 0
		let a1, b1, c1, a2, b2;
		if (edgeDeltaX === 0) {
			// x = -c1
			a1 = 1;
			b1 = 0;
			c1 = -edgeX1;
			// y = -c2
			a2 = 0;
			b2 = 1;
		} else if (edgeDeltaY === 0) {
			// y = -c1
			a1 = 0;
			b1 = 1;
			c1 = -edgeY1;
			// x = -c2
			a2 = 1;
			b2 = 0;
		} else {
			a1 = -edgeDeltaY / edgeDeltaX;
			b1 = 1;
			c1 = -a1 * edgeX1 - edgeY1;
			a2 = edgeDeltaX / edgeDeltaY;
			b2 = 1;
		}

		let maxDistanceSq = 0;
		let x1, y1, x2, y2;
		const candidateX = [0, edgeX1, edgeX2];
		const candidateY = [0, edgeY1, edgeY2];
		for (let i = 0; i < numPoints; i++) {
			const targetX = pointsX[i];
			const targetY = pointsY[i];
			const c2 = -a2 * targetX - b2 * targetY;
			const denominator = a1 * b2 - a2 * b1;
			const intersectX = (b1 * c2 - b2 * c1) / denominator;
			const intersectY = (c1 * a2 - c2 * a1) / denominator;
			candidateX[0] = intersectX;
			candidateY[0] = intersectY;
			let candidate = 0;

			if (edgeX1 < edgeX2) {
				if (intersectX > edgeX2) {
					candidate = 2;
				} else if (intersectX < edgeX1) {
					candidate = 1;
				}
			} else if (edgeX2 < edgeX1) {
				if (intersectX > edgeX1) {
					candidate = 1;
				} else if (intersectX < edgeX2) {
					candidate = 2;
				}
			} else {
				if (edgeY1 < edgeY2) {
					if (intersectY > edgeY2) {
						candidate = 2;
					} else if (intersectY < edgeY1) {
						candidate = 1;
					}
				} else if (edgeY2 < edgeY1) {
					if (intersectY > edgeY1) {
						candidate = 1;
					} else if (intersectY < edgeY2) {
						candidate = 2;
					}
				}
			}

			const measureX = candidateX[candidate];
			const measureY = candidateY[candidate];
			const gradientDeltaX = targetX - measureX;
			const gradientDeltaY = targetY - measureY;
			const distanceSq = gradientDeltaX * gradientDeltaX + gradientDeltaY * gradientDeltaY;
			if (distanceSq > maxDistanceSq) {
				x1 = intersectX;
				y1 = intersectY;
				x2 = targetX;
				y2 = targetY;
				maxDistanceSq = distanceSq;
			}
		}

		const colour1 = this.colourMorph1.interpolate(morph, interpolation);
		const colour2 = this.colourMorph2.interpolate(morph, interpolation);
		const gradient = context.createLinearGradient(x1, y1, x2, y2);
		gradient.addColorStop(0, colour1);
		gradient.addColorStop(1, colour2);
		return gradient;
	}

}

class StrokeStyle {
	constructor() {
		this.width = undefined;
		this.dash = undefined;
		this.dashOffset = undefined;
		this.start = 0;
		this.end = 1;
	}
}

class StrokeMorph {

	constructor() {
		this.colourMorph = undefined;
		this.startStyle = new StrokeStyle();
		this.endStyle = new StrokeStyle();
	}

	setColour(colourMorph) {
		this.colourMorph = colourMorph;
	}

	setWidth(start, end = start) {
		this.startStyle.width = start;
		this.endStyle.width = end;
	}

	setDash(start, end = start) {
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

class Morph {

	constructor(
		polygon1, polygon2, fillMorph, strokeMorph, blendMode = 'source-over', maxRotation = DEFAULT_MAX_ROTATION
	) {
		polygon2.resize(polygon1.size);
		polygon2.rotate(polygon1, maxRotation);

		this.polygon1 = polygon1;
		this.polygon2 = polygon2;
		this.numFrames = 64;
		this.interpolationStep = 2 ** -6;
		this.translateX = 0;
		this.translateY = 0;
		this.scale = 1;
		this.rotation = 0;
		this.pointsX = undefined;
		this.pointsY = undefined;
		this.fillMorph = fillMorph;
		this.fillStyle = 'black';
		this.strokeMorph = strokeMorph;
		this.stroke = undefined;
		this.blendMode = blendMode;
	}

	setSpeed(speed) {
		// Technically the actual number of frames is one more than this.
		const numFrames = Math.ceil(moveLength(this.polygon1, this.polygon2) / speed);
		this.numFrames = numFrames;
		this.interpolationStep = 1 / numFrames;
	}

	interpolate(interpolation, context) {
		const polygon1 = this.polygon1;
		const polygon2 = this.polygon2;
		this.translateX = polygon2.centreX * interpolation + polygon1.centreX * (1 - interpolation);
		this.translateY = polygon2.centreY * interpolation + polygon1.centreY * (1 - interpolation);
		this.scale = polygon2.size / polygon1.size * interpolation + 1 - interpolation;
		this.rotation = -polygon2.rotation * interpolation;

		const numPoints = polygon2.numPoints;
		const pointsX = new Array(numPoints);
		const pointsY = new Array(numPoints);
		for (let i = 0; i < numPoints; i++) {
			pointsX[i] = polygon2.rotatedX[i] - (1 - interpolation) * polygon2.offsetsX[i];
			pointsY[i] = polygon2.rotatedY[i] - (1 - interpolation) * polygon2.offsetsY[i];
		}
		this.pointsX = pointsX;
		this.pointsY = pointsY;

		if (this.fillMorph) {
			this.fillStyle = this.fillMorph.interpolate(this, interpolation, context);
		}

	}

	transform(context) {
		context.translate(this.translateX, this.translateY);
		context.scale(this.scale, this.scale);
		context.rotate(this.rotation);
		context.fillStyle = this.fillStyle;
	}

	animate(context, mainCallback, beforeCallback = undefined) {
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
			me.interpolate(interpolation, context);

			const canvas = context.canvas;
			context.clearRect(0, 0, canvas.width, canvas.height);
			if (beforeCallback) {
				beforeCallback(context, me, interpolation);
			}
			context.save();
			me.transform(context);
			context.globalCompositeOperation = me.blendMode;
			mainCallback(context, me, interpolation);
			context.restore();

			if (interpolation < 1) {
				requestAnimationFrame(drawFrame);
			}
		}

		requestAnimationFrame(drawFrame);
	}

	overlay(context, callback) {
		this.interpolate(0, context);
		context.save();
		this.transform(context);
		callback(context, this, 0);
		context.restore();

		for (
			let interpolation = this.interpolationStep;
			interpolation <= 1;
			interpolation += this.interpolationStep
		) {
			this.interpolate(interpolation, context);
			context.save();
			this.transform(context);
			context.globalCompositeOperation = this.blendMode;
			callback(context, this, interpolation);
			context.restore();
		}
	}

}


function linearInterpolation(startX, startY, endX, endY) {
	const xDistance = endX - startX;
	const yDistance = endY - startY;
	function interpolate(t) {
		const x = startX + xDistance * t;
		const y = startY + yDistance * t;
		return [x, y];
	}
	return interpolate;
}

class Path2D {

	static generate(pathGenerator, startT = 0, endT = 1) {
		let [startX, startY] = pathGenerator(startT);
		startX = Math.round(startX);
		startY = Math.round(startY);

		let [endX, endY] = pathGenerator(endT);
		endX = Math.round(endX);
		endY = Math.round(endY);
		let path;
		if (startX !== endX || startY !== endY) {
			// open curve
			path = Path2D.#interpolateSegment(
				pathGenerator, startX, startY, startT, endX, endY, endT
			);
			path.prependPoint(startX, startY, startT);
			path.appendPoint(endX, endY, endT);
		} else {
			// closed curve
			let lowMidT = endT;
			let highMidT = startT;
			let lowMidX, lowMidY, highMidX, highMidY;
			do {
				lowMidT = 0.5 * (startT + lowMidT);
				highMidT = 0.5 * (highMidT + endT);

				[lowMidX, lowMidY] = pathGenerator(lowMidT);
				lowMidX = Math.round(lowMidX);
				lowMidY = Math.round(lowMidY);

				[highMidX, highMidY] = pathGenerator(highMidT);
				highMidX = Math.round(highMidX);
				highMidY = Math.round(highMidY);
			} while (lowMidX === startX && lowMidY === startY && highMidX === endX && highMidY === endY);
			path = Path2D.#interpolateSegment(
				pathGenerator, startX, startY, startT, lowMidX, lowMidY, lowMidT
			);
			path.prependPoint(startX, startY, startT);
			path.appendPoint(lowMidX, lowMidY, lowMidT);
			let path2;
			if (lowMidT !== highMidT) {
				path2 = Path2D.#interpolateSegment(
					pathGenerator, lowMidX, lowMidY, lowMidT, highMidX, highMidY, highMidT
				);
				path.appendPath(path2);
				path.appendPoint(highMidX, highMidY, highMidT);
			}
			path2 = Path2D.#interpolateSegment(
				pathGenerator, highMidX, highMidY, highMidT, endX, endY, endT
			);
			path.appendPath(path2);
			path.appendPoint(endX, endY, endT);
		}
		return path;
	}

	static #interpolateSegment(pathGenerator, startX, startY, startT, endX, endY, endT) {
		const xDistance = Math.abs(endX - startX);
		const yDistance = Math.abs(endY - startY);

		if (xDistance <= 1 && yDistance <= 1) {
			// No intermediate points needed. The start and end points of the current segment
			// are enough.
			return new Path2D([], [], []);
		}

		const midT = 0.5 * (startT + endT);
		let [midX, midY] = pathGenerator(midT);
		midX = Math.round(midX);
		midY = Math.round(midY);

		const path = Path2D.#interpolateSegment(
			pathGenerator, startX, startY, startT, midX, midY, midT
		);
		path.appendPoint(midX, midY, midT);
		const path2 = Path2D.#interpolateSegment(
			pathGenerator, midX, midY, midT, endX, endY, endT
		);
		path.appendPath(path2);
		return path;
	}

	constructor(xValues, yValues, tValues) {
		this.xValues = xValues;
		this.yValues = yValues;
		this.tValues = tValues;
	}

	get length() {
		return this.tValues.length;
	}

	getPoint(searchT) {
		const tValues = this.tValues;
		if (searchT === tValues[0]) {
			return [ this.xValues[0], this.yValues[0] ];
		}
		let upperIndex = tValues.length - 1;
		if (searchT === tValues[upperIndex]) {
			return [ this.xValues[upperIndex], this.yValues[upperIndex] ];
		}

		let lowerIndex = 0;
		let resultIndex;
		while (true) {
			const midIndex = (lowerIndex + upperIndex) >> 1;
			const midValue = tValues[midIndex];
			if (searchT === midValue) {
				resultIndex = midIndex;
				break;
			} else if (searchT < midValue) {
				if (lowerIndex = midIndex - 1) {
					if (searchT - tValues[lowerIndex] <= midValue - searchT) {
						resultIndex = lowerIndex;
					} else {
						resultIndex = midIndex;
					}
					break;
				}
				upperIndex = midIndex - 1;
			} else {
				// searchT > midValue
				if (upperIndex = midIndex + 1) {
					if (searchT - midValue <= tValues[upperIndex] - searchT) {
						resultIndex = midIndex;
					} else {
						resultIndex = upperIndex;
					}
					break;
				}
				lowerIndex = midIndex + 1;
			}
		} // end loop
		return [ this.xValues[resultIndex], this.yValues[resultIndex] ];
	}

	prependPoint(x, y, t) {
		this.xValues.unshift(x);
		this.yValues.unshift(y);
		this.tValues.unshift(t);
	}

	appendPoint(x, y, t) {
		this.xValues.push(x);
		this.yValues.push(y);
		this.tValues.push(t);
	}

	appendPath(path) {
		const numPoints = this.tValues.length;
		this.xValues.splice(numPoints, 0, ...path.xValues);
		this.yValues.splice(numPoints, 0, ...path.yValues);
		this.tValues.splice(numPoints, 0, ...path.tValues);
	}

	/**Computes some statistic about the rate of change in the curve and returns that statistic.
	 *
	 * Examples:
	 *
	 * initialValue	func			Result
	 * ---------------------------------
	 * Infinity			Math.min		Value of delta t to use to animate the slowest part of the
	 * 									curve at 1 pixel per frame.
	 * 0					Math.max		Value of delta t to use to animate the fastest part of the
	 * 									curve at 1 pixel per frame.
	 */
	deltaTCalculation(initialValue, func) {
		const xValues = this.xValues;
		const yValues = this.yValues;
		const tValues = this.tValues;
		const numPoints = tValues.length;
		let prevX = xValues[0];
		let prevY = yValues[0];
		let prevT = tValues[0];
		let deltaTStatistic = initialValue;
		for (let i = 1; i < numPoints; i++) {
			const x = xValues[i];
			const y = yValues[i];
			const t = tValues[i];
			const xDistance = Math.abs(x - prevX);
			const yDistance = Math.abs(y - prevY);
			const deltaSpace = Math.max(xDistance, yDistance);
			if (deltaSpace > 0) {
				const deltaT = (t - prevT) / deltaSpace;
				deltaTStatistic = func(deltaTStatistic, deltaT);
				prevX = x;
				prevY = y;
				prevT = t;
			}
		}
		return deltaTStatistic;
	}

}


function drawSourceShapes(context, morph, interpolation) {
	const polygon1 = morph.polygon1;
	const polygon2 = morph.polygon2;
	context.globalAlpha = 0.39;

	polygonPath(context, polygon1.pointsX, polygon1.pointsY);
	context.fillStyle = 'red';
	context.fill();

	polygonPath(context, polygon2.pointsX, polygon2.pointsY);
	context.fillStyle = '#ccc';
	context.fill();
}

function drawInterpolatedShape(context, morph, interpolation) {
	polygonPath(context, morph.pointsX, morph.pointsY);
	context.fill();
}

function drawFaded(context, morph, interpolation) {
	const alpha = Math.max(Math.round(0.9 * 255 / morph.numFrames), 1) / 255;
	context.globalAlpha = alpha;

	polygonPath(context, morph.pointsX, morph.pointsY);
	context.fill();
}


const canvas = document.getElementById('canvas');
const context = canvas.getContext('2d');
canvas.width = window.innerWidth;
canvas.height = window.innerHeight;

const Mode = Object.freeze({
	ANIMATE: 0,
	OVERLAY: 1,
});

const parameters = new URLSearchParams(document.location.search);
let speed = parseFloat(parameters.get('speed'));
let mode, fillMorph;

const fillStr = parameters.get('fill');
if (fillStr) {
	// Use %25 to write % inside a URL and %2B to write +.
	fillMorph = new ColourMorph(fillStr);
}

switch (parameters.get('render')) {
case 'overlay':
	mode = Mode.OVERLAY;
	speed ||= 50;
	if (!fillMorph) {
		fillMorph = new ColourMorph('rgb(calc(255-510*t), min(255*t, 255-255*t), calc(510*t-255))');
	}
	break;
default:
	mode = Mode.ANIMATE;
	speed ||= 6;
	if (!fillMorph) {
		fillMorph = new ConstantColour('lime');
	}
}

const blendMode = parameters.get('blend') || 'source-over';

const numVertices = parseInt(parameters.get('vertices')) || 5;
let polygon1 = randomPolygon(numVertices);
let polygon2 = randomPolygon(numVertices);
if (mode === Mode.OVERLAY && polygon2.size > polygon1.size) {
	const temp = polygon1;
	polygon1 = polygon2;
	polygon2 = temp;
}

const fillStr2 = parameters.get('fill2');
if (fillStr2) {
	const fromVertex = 0;
	const toVertex = Math.trunc(numVertices / 2);
	const toColour = new ColourMorph(fillStr2);
	//fillMorph = new VertexGradientMorph(fillMorph, toColour, fromVertex, toVertex);
	fillMorph = new EdgeGradientMorph(fillMorph, toColour, fromVertex, toVertex);
}

let maxRotation = parseInt(parameters.get('max_rotation'));
if (Number.isFinite(maxRotation)) {
	maxRotation *= 180 / Math.PI;
} else {
	maxRotation = DEFAULT_MAX_ROTATION;
}

const morph = new Morph(polygon1, polygon2, fillMorph, undefined, blendMode, maxRotation);
morph.setSpeed(speed);

switch (mode) {
case Mode.ANIMATE:
	// Draw as an animation.
	morph.animate(context, drawInterpolatedShape, drawSourceShapes);
	break;
case Mode.OVERLAY:
	// Draw with successive frames overlaid on top of each other.
	morph.overlay(context, drawFaded);
}

