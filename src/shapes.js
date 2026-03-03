// Code for modelling shapes.

import {Geometry, Sort} from './math.js';

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
		[pointsX, pointsY] = Geometry.orderPolygonPoints(pointsX, pointsY);
		const [centreX, centreY] = Geometry.centroid(pointsX, pointsY);

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

		this.pointsX = pointsX;
		this.pointsY = pointsY;
		this.centreX = centreX;
		this.centreY = centreY;
		this.deltaX = deltaX;
		this.deltaY = deltaY;
		this.angles = angles;
		this.radii = radii;
		this.size = Math.sqrt(area / Math.PI);
		this.offsetsX = new Array(numPoints);
		this.offsetsY = new Array(numPoints);
		this.resetTransform();
	}

	get numPoints() {
		return this.pointsX.length;
	}

	resetTransform() {
		const deltaX = this.deltaX;
		const deltaY = this.deltaY;
		this.resizedX = deltaX.slice();
		this.resizedY = deltaY.slice();
		this.scaleFactor = 1;
		this.rotatedX = deltaX.slice();
		this.rotatedY = deltaY.slice();
		this.rotation = 0;
		// The offsets that remain after translation, resizing and rotating.
		this.offsetsX.fill(0);
		this.offsetsY.fill(0);
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
		this.scaleFactor = scale;
	}

	/**Adds a new vertex without recalculating the centre point or average radius.
	 * @param {number} index The index of the vertex which should end up immediately preceding
	 * the new vertex.
	 */
	addVertex(index, rotatedX, rotatedY) {
		this.rotatedX.splice(index, 0, rotatedX);
		this.rotatedY.splice(index, 0, rotatedY);
		this.offsetsX.splice(index, 0, 0);
		this.offsetsY.splice(index, 0, 0);

		const cos = Math.cos(-this.rotation);
		const sin = Math.sin(-this.rotation);
		const resizedX = rotatedX * cos - rotatedY * sin;
		const resizedY = rotatedX * sin + rotatedY * cos;
		this.resizedX.splice(index, 0, resizedX);
		this.resizedY.splice(index, 0, resizedY);

		const angle = Math.atan2(resizedY, resizedX);
		this.angles.splice(index, 0, angle);
		const radius = (resizedX / this.scaleFactor) / Math.cos(angle);
		this.radii.splice(index, 0, radius);
		const deltaX = resizedX / this.scaleFactor;
		const deltaY = resizedY / this.scaleFactor;
		this.deltaX.splice(index, 0, deltaX);
		this.deltaY.splice(index, 0, deltaY);
		this.pointsX.splice(index, 0, deltaX + this.centreX);
		this.pointsY.splice(index, 0, deltaY + this.centreY);

	}

	/**Computes a new set of points by rotating the scaled shape to align with the target shape
	 * as best as possible and stores the rectangular coordinates of the new points in
	 * this.rotatedX and this.rotatedY, along with the residual displacement vectors in
	 * this.offsetsX and this.offsetsY.
	 */
	rotate(shape2, maxRotation) {
		const numPoints = this.numPoints;
		const numPoints2 = shape2.numPoints;
		const minNumPoints = Math.min(numPoints, numPoints2);

		let selectedRotation = 0;
		let minError = Infinity;
		let minErrorAlignmentIndex = 0;
		let minErrorChoice, minErrorChoice2;

		for (let i = 0; i < numPoints2; i++) {
			for (let choice of Geometry.choices(numPoints, minNumPoints)) {
				let rotation = shape2.angles[i] - this.angles[choice[0]];
				if (rotation < -Math.PI) {
					rotation += 2 * Math.PI;
				} else if (rotation > Math.PI) {
					rotation -= 2 * Math.PI;
				}
				const performRotation = Math.abs(rotation) <= maxRotation;
				const cos = Math.cos(rotation);
				const sin = Math.sin(rotation);

				for (let choice2 of Geometry.choices(numPoints2, minNumPoints)) {
					const alignmentIndex = choice2.indexOf(i);
					if (alignmentIndex === -1) {
						continue;
					}
					let sumOfSquares = 0;
					for (let j = 0; j < minNumPoints; j++) {
						const sourceIndex = choice[j];
						const x = this.resizedX[sourceIndex];
						const y = this.resizedY[sourceIndex];
						let transformedX = x, transformedY = y;
						if (performRotation) {
							transformedX = x * cos - y * sin;
							transformedY = x * sin + y * cos;
						}
						const targetIndex = choice2[(alignmentIndex + j) % minNumPoints];
						const targetX = shape2.resizedX[targetIndex];
						const targetY = shape2.resizedY[targetIndex];
						const distanceX = targetX - transformedX;
						const distanceY = targetY - transformedY;
						sumOfSquares += distanceX * distanceX + distanceY * distanceY;
					}
					if (sumOfSquares < minError) {
						minError = sumOfSquares;
						minErrorChoice = choice;
						minErrorChoice2 = choice2;
						minErrorAlignmentIndex = alignmentIndex;
						selectedRotation = rotation;
					}
				}	// end for each choice of vertices in Shape 2
			}		// end for each choice of vertices in Shape 1
		}			// end for each vertex of Shape 2 to align to
		this.rotation = selectedRotation;
		this.rotatedX = new Array(numPoints);
		this.rotatedY = new Array(numPoints);
		const cos = Math.cos(selectedRotation);
		const sin = Math.sin(selectedRotation);
		for (let i = 0; i < numPoints; i++) {
			this.rotatedX[i] = this.resizedX[i] * cos - this.resizedY[i] * sin;
			this.rotatedY[i] = this.resizedX[i] * sin + this.resizedY[i] * cos;
		}
		for (let i = 0; i < minNumPoints; i++) {
			const source = minErrorChoice[i];
			const target = minErrorChoice2[(minErrorAlignmentIndex + i) % minNumPoints];
			this.offsetsX[source] = this.rotatedX[source] - shape2.resizedX[target];
			this.offsetsY[source] = this.rotatedY[source] - shape2.resizedY[target];
		}

		let numInserted = 0;
		if (numPoints < numPoints2) {
			let lastTarget = minErrorChoice2[minErrorAlignmentIndex];
			let edgeX1 = this.rotatedX[0];
			let edgeY1 = this.rotatedY[0];
			for (let i = 1; i <= minNumPoints; i++) {
				let source2 = (i + numInserted) % (numPoints + numInserted);
				const edgeX2 = this.rotatedX[source2];
				const edgeY2 = this.rotatedY[source2];
				const [edgeDeltaX, a1, b1, c1, a2, b2] = Geometry.projectionOntoLine(
					edgeX1, edgeY1, edgeX2, edgeY2
				);
				const nextTarget = minErrorChoice2[(minErrorAlignmentIndex + i) % minNumPoints];
				let numIntermediate = nextTarget - lastTarget - 1;
				if (numIntermediate < 0) {
					numIntermediate += numPoints2;
				}
				for (let j = 1; j <= numIntermediate; j++) {
					const target = (lastTarget + j) % numPoints2;
					const targetX = shape2.resizedX[target];
					const targetY = shape2.resizedY[target];
					let [projectX, projectY] = Geometry.projectOntoLine(targetX, targetY, a1, b1, c1, a2, b2);
					[projectX, projectY] = Geometry.constrainToLineSegment(
						projectX, projectY, edgeX1, edgeY1, edgeX2, edgeY2
					);
					this.addVertex(source2, projectX, projectY);
					this.offsetsX[source2] = projectX - targetX;
					this.offsetsY[source2] = projectY - targetY;
					source2++;
					numInserted++;
				}
				lastTarget = nextTarget;
				edgeX1 = edgeX2;
				edgeY1 = edgeY2;
			}

		} else if (numPoints > numPoints2) {

			let lastSource = minErrorChoice[0];
			let edgeX1 = shape2.resizedX[minErrorAlignmentIndex];
			let edgeY1 = shape2.resizedY[minErrorAlignmentIndex];
			for (let i = 1; i <= minNumPoints; i++) {
				let target2 = (minErrorAlignmentIndex + i + numInserted) %
					(numPoints2 + numInserted);
				const edgeX2 = shape2.resizedX[target2];
				const edgeY2 = shape2.resizedY[target2];
				const [edgeDeltaX, a1, b1, c1, a2, b2] = Geometry.projectionOntoLine(
					edgeX1, edgeY1, edgeX2, edgeY2
				);
				const nextSource = minErrorChoice[i % minNumPoints];
				let numIntermediate = nextSource - lastSource - 1;
				if (numIntermediate < 0) {
					numIntermediate += numPoints;
				}
				for (let j = 1; j <= numIntermediate; j++) {
					const source = (lastSource + j) % numPoints;
					const sourceX = this.rotatedX[source];
					const sourceY = this.rotatedY[source];
					let [projectX, projectY] = Geometry.projectOntoLine(sourceX, sourceY, a1, b1, c1, a2, b2);
					[projectX, projectY] = Geometry.constrainToLineSegment(
						projectX, projectY, edgeX1, edgeY1, edgeX2, edgeY2
					);
					shape2.addVertex(target2, projectX, projectY);
					this.offsetsX[source] = sourceX - projectX;
					this.offsetsY[source] = sourceY - projectY;
					if (target2 <= minErrorAlignmentIndex) {
						minErrorAlignmentIndex++;
					}
					target2++;
					numInserted++;
				}
				lastSource = nextSource;
				edgeX1 = edgeX2;
				edgeY1 = edgeY2;
			}

		}

	}

}

function randomPolygonPoints(numPoints, maxWidth, maxHeight) {
	const pointsX = new Array(numPoints);
	const pointsY = new Array(numPoints);
	for (let i = 0; i < numPoints; i++) {
		pointsX[i] = Math.trunc(Math.random() * maxWidth);
		pointsY[i] = Math.trunc(Math.random() * maxHeight);
	}
	return [pointsX, pointsY];
}

/**Returns a value from a probability distribution which has a 50% chance of being between
 * minFraction and 1 and a 50% of being between 1 and maxValue.
*/
function randomFraction(minFraction, maxValue) {
	const r = Math.random() * 2 * (maxValue - 1);
	if (r >= maxValue - 1) {
		return r - (maxValue - 1) + 1;
	}
	return r / (maxValue - 1) * (1 - minFraction) + minFraction;
}

/**Creates a random cyclic polygon centred on (0,0) with radius 1 and order 2 reflective
 * symmetry.
 */
function cyclicPolygonPoints(
	numPoints, radiusX, radiusY = radiusX, offsetX = 0, offsetY = 0, minAngleProportion = 0.5,
	maxAngleProportion = 1.5, minRotation = 0, maxRotation = 2 * Math.PI
) {
	const angles = new Array(numPoints);
	angles[0] = 0;
	const regularAngle = 2 * Math.PI / numPoints;
	const numFreeAngles = Math.ceil(0.5 * numPoints) - 1;
	let amountLeft = Math.PI;
	if ((numPoints & 1) === 0) {
		angles[0.5 * numPoints] = Math.PI;
		for (let i = 1; i <= numFreeAngles; i++) {
			const equalAmount = amountLeft / (numFreeAngles - i + 2);
			const minAmount = Math.max(
				amountLeft - (numFreeAngles - i + 1) * maxAngleProportion * regularAngle,
				minAngleProportion * regularAngle
			);
			const maxAmount = Math.min(
				amountLeft - (numFreeAngles - i + 1) * minAngleProportion * regularAngle,
				maxAngleProportion * regularAngle
			);
			const amount = Math.min(Math.max(
				randomFraction(minAngleProportion, maxAngleProportion) * equalAmount,
				minAmount), maxAmount);
			amountLeft -= amount;
			angles[i] = Math.PI - amountLeft;
		}
	} else {
		for (let i = 1; i <= numFreeAngles; i++) {
			const equalAmount = amountLeft / (numFreeAngles - i + 1.5);
			const minAmount = Math.max(
				amountLeft - (numFreeAngles - i + 0.5) * maxAngleProportion * regularAngle,
				minAngleProportion * regularAngle
			);
			const maxAmount = Math.min(
				amountLeft - (numFreeAngles - i + 0.5) * minAngleProportion * regularAngle,
				maxAngleProportion * regularAngle
			);
			const amount = Math.min(Math.max(
				randomFraction(minAngleProportion, maxAngleProportion) * equalAmount,
				minAmount), maxAmount);
			amountLeft -= amount;
			angles[i] = Math.PI - amountLeft;
		}
	}
	for (let i = 0; i < numFreeAngles; i++) {
		angles[numPoints - 1 - i] = 2 * Math.PI - angles[i + 1];
	}
	const rotation = Math.random() * (maxRotation - minRotation) + minRotation;
	const pointsX = new Array(numPoints);
	const pointsY = new Array(numPoints);
	for (let i = 0; i < numPoints; i++) {
		const angle = angles[i] + rotation;
		pointsX[i] = Math.round(radiusX * Math.cos(angle) + offsetX);
		pointsY[i] = Math.round(radiusY * Math.sin(angle) + offsetY);
	}
	return [pointsX, pointsY]
}

class Path {

	static generate(
		pathGenerator, startT = 0, endT = 1, parameter = 0, roundX = 1, roundY = roundX
	) {
		let [startX, startY] = pathGenerator(startT, parameter);
		startX = Math.round(startX / roundX) * roundX;
		startY = Math.round(startY / roundY) * roundY;

		let [endX, endY] = pathGenerator(endT, parameter);
		endX = Math.round(endX / roundX) * roundX;
		endY = Math.round(endY / roundY) * roundY;

		if (startX !== endX || startY !== endY) {
			// Open curve
			const path = Path.#interpolateSegment(
				pathGenerator, parameter, startX, startY, startT, endX, endY, endT
			);
			path.#prependPoint(startX, startY, startT);
			path.#appendPoint(endX, endY, endT);
			return path;
		}

		// Closed curve
		let lowMidT = 0.5 * (startT + endT);
		let highMidT = lowMidT;

		let [lowMidX, lowMidY] = pathGenerator(lowMidT, parameter);
		lowMidX = Math.round(lowMidX / roundX) * roundX;
		lowMidY = Math.round(lowMidY / roundY) * roundY;
		let highMidX = lowMidX;
		let highMidY = lowMidY;

		const MAX_ITERATIONS = 6;	// Safety measure
		let i = 0;
		while (lowMidX === startX && lowMidY === startY && i < MAX_ITERATIONS) {
			lowMidT = 0.5 * (startT + lowMidT);
			[lowMidX, lowMidY] = pathGenerator(lowMidT, parameter);
			lowMidX = Math.round(lowMidX / roundX) * roundX;
			lowMidY = Math.round(lowMidY / roundY) * roundY;
			i++;
		}
		i = 0;
		while (highMidX === endX && highMidY === endY && i < MAX_ITERATIONS) {
			highMidT = 0.5 * (highMidT + endT);
			[highMidX, highMidY] = pathGenerator(highMidT, parameter);
			highMidX = Math.round(highMidX / roundX) * roundX;
			highMidY = Math.round(highMidY / roundY) * roundY;
			i++;
		}

		// Range [startT, lowMidT]
		let path = Path.#interpolateSegment(
			pathGenerator, parameter, startX, startY, startT, lowMidX, lowMidY, lowMidT, roundX,
			roundY
		);
		path.#prependPoint(startX, startY, startT);
		path.#appendPoint(lowMidX, lowMidY, lowMidT);

		let path2;

		// Range (lowMidT, highMidT]
		if (lowMidT !== highMidT) {
			path2 = Path.#interpolateSegment(
				pathGenerator, parameter, lowMidX, lowMidY, lowMidT, highMidX, highMidY, highMidT,
				roundX, roundY
			);
			path.#appendPath(path2);
			path.#appendPoint(highMidX, highMidY, highMidT);
		}

		// Range (highMidT, endT]
		path2 = Path.#interpolateSegment(
			pathGenerator, parameter, highMidX, highMidY, highMidT, endX, endY, endT, roundX,
			roundY
		);
		path.#appendPath(path2);
		path.#appendPoint(endX, endY, endT);
		return path;
	}

	static #interpolateSegment(
		pathGenerator, parameter, startX, startY, startT, endX, endY, endT, roundX, roundY
	) {
		const xDistance = Math.abs(endX - startX);
		const yDistance = Math.abs(endY - startY);

		if (xDistance <= roundX && yDistance <= roundY) {
			// No intermediate points needed. The start and end points of the current segment
			// are enough.
			return new Path([], [], []);
		}

		const midT = 0.5 * (startT + endT);
		let [midX, midY] = pathGenerator(midT, parameter);
		midX = Math.round(midX / roundX) * roundX;
		midY = Math.round(midY / roundY) * roundY;

		const path = Path.#interpolateSegment(
			pathGenerator, parameter, startX, startY, startT, midX, midY, midT, roundX, roundY
		);
		path.#appendPoint(midX, midY, midT);
		const path2 = Path.#interpolateSegment(
			pathGenerator, parameter, midX, midY, midT, endX, endY, endT, roundX, roundY
		);
		path.#appendPath(path2);
		return path;
	}

	constructor(xValues, yValues, tValues) {
		this.pointsX = xValues;
		this.pointsY = yValues;
		this.pointsT = tValues;
	}

	get numPoints() {
		return this.pointsT.length;
	}

	getIndex(searchT) {
		let lowerIndex = 0;
		let upperIndex = this.pointsT.length - 1;

		while (true) {
			const midIndex = (lowerIndex + upperIndex) >> 1;
			const midT = this.pointsT[midIndex];

			if (searchT === midT) {
				return midIndex;
			}

			if (searchT < midT) {

				if (lowerIndex === midIndex - 1) {
					if (searchT - this.pointsT[lowerIndex] < midT - searchT) {
						return lowerIndex;
					} else {
						return midIndex;
					}
				}
				upperIndex = midIndex - 1;

			} else {

				// searchT > midT
				if (upperIndex === midIndex + 1) {
					if (searchT - midT < this.pointsT[upperIndex] - searchT) {
						return midIndex;
					} else {
						return upperIndex;
					}
				}
				lowerIndex = midIndex + 1;
			}

		} // end loop
	}

	#prependPoint(x, y, t) {
		this.pointsX.unshift(x);
		this.pointsY.unshift(y);
		this.pointsT.unshift(t);
	}

	#appendPoint(x, y, t) {
		this.pointsX.push(x);
		this.pointsY.push(y);
		this.pointsT.push(t);
	}

	#appendPath(path) {
		const numPoints = this.pointsT.length;
		this.pointsX.splice(numPoints, 0, ...path.pointsX);
		this.pointsY.splice(numPoints, 0, ...path.pointsY);
		this.pointsT.splice(numPoints, 0, ...path.pointsT);
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
	statistic(initialValue, func) {
		let prevT = this.pointsT[0];
		let statistic = initialValue;
		for (let i = 1; i < this.numPoints; i++) {
			const t = this.pointsT[i];
			statistic = func(statistic, t - prevT);
			prevT = t;
		}
		return statistic;
	}

}

const RandomShape = {
	polygonPoints: randomPolygonPoints,
	cyclicPolygonPoints: cyclicPolygonPoints,
}

export {
	Shape,
	RandomShape,
	Path,
};
