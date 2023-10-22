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
	for (let i = 0; i < numPoints; i++) {
		maxTranslation = Math.max(
			maxTranslation,
			Math.abs(translateX + polygon2.offsetsX[i]),
			Math.abs(translateY + polygon2.offsetsY[i])
		);

		if (absRotation > 0) {
			const radius = polygon2.radii[i];
			const scaledRadius = radius * inverseScale;
			// Arc length of an Archimedean spiral.
			const radialGrowth = (scaledRadius - radius) / absRotation;
			let squareRoot = Math.sqrt(scaledRadius * scaledRadius + radialGrowth * radialGrowth);
			const firstTerm =  0.5 * squareRoot * scaledRadius / radialGrowth;
			const secondTerm = 0.5 * radialGrowth * Math.log(squareRoot + scaledRadius);
			squareRoot = Math.sqrt(radius * radius + radialGrowth * radialGrowth);
			const thirdTerm = 0.5 * squareRoot * radius / radialGrowth;
			const fourthTerm = 0.5 * radialGrowth * Math.log(squareRoot + radius);
			const arcLength = (firstTerm + secondTerm) - (thirdTerm + fourthTerm);
			maxArcLength = Math.max(maxArcLength, arcLength);
		}
	}
	return Math.max(maxTranslation, maxArcLength);
}

class Morph {

	constructor(polygon1, polygon2, maxRotation = DEFAULT_MAX_ROTATION) {
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
	}

	setSpeed(speed) {
		// Technically the actual number of frames is one more than this.
		const numFrames = Math.ceil(moveLength(this.polygon1, this.polygon2) / speed);
		this.numFrames = numFrames;
		this.interpolationStep = 1 / numFrames;
	}

	interpolate(interpolation) {
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
	}

	transform(context) {
		context.translate(this.translateX, this.translateY);
		context.scale(this.scale, this.scale);
		context.rotate(this.rotation);
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
			me.interpolate(interpolation);

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

	overlay(context, callback) {
		for (let interpolation = 0; interpolation <= 1; interpolation += this.interpolationStep) {
			this.interpolate(interpolation);
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
	context.globalAlpha = 0.39;

	polygonPath(context, polygon1.pointsX, polygon1.pointsY);
	context.fillStyle = 'red';
	context.fill();

	polygonPath(context, polygon2.pointsX, polygon2.pointsY);
	context.fillStyle = 'blue';
	context.fill();
}

function drawInterpolatedShape(context, morph, interpolation) {
	polygonPath(context, morph.pointsX, morph.pointsY);
	context.fillStyle = 'lime';
	context.fill();
}

function drawFaded(context, morph, interpolation) {
	context.globalAlpha = Math.max(
		Math.round(0.9 * 255 / morph.numFrames),
		1
	) / 255;
	let startComponent = interpolation < 0.5 ? Math.round((1 - 2 * interpolation) * 255) : 0;
	let endComponent = interpolation > 0.5 ? Math.round(2 * (interpolation - 0.5) * 255) : 0;
	let midComponent = Math.round(Math.min(interpolation, 1 - interpolation) * 255);
	context.fillStyle = `rgb(${startComponent}, ${midComponent}, ${endComponent})`;

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
let mode;

switch (parameters.get('render')) {
case 'overlay':
	mode = Mode.OVERLAY;
	speed ||= 50;
	break;
default:
	mode = Mode.ANIMATE;
	speed ||= 6;
}

let numVertices = parseInt(parameters.get('vertices')) || 5;


let polygon1 = randomPolygon(numVertices);
let polygon2 = randomPolygon(numVertices);
if (mode === Mode.OVERLAY && polygon2.size > polygon1.size) {
	const temp = polygon1;
	polygon1 = polygon2;
	polygon2 = temp;
}

let maxRotation = parseInt(parameters.get('max_rotation'));
if (maxRotation === null) {
	maxRotation = DEFAULT_MAX_ROTATION;
} else {
	maxRotation *= 180 / Math.PI;
}

const morph = new Morph(polygon1, polygon2, maxRotation);
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

