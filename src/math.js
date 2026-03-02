/* Functions for finding general mathematical qualities required to compute shape morphs, with
 * a focus on geometry.
 */

//Sorting functions

const Sort = {
	numericAscending: (a, b) => a - b,
	numericDescending: (a, b) => b - a,
}


// General maths, e.g. combinatorics

function gcd(a, b) {
	while (b !== 0) {
		const temp = b;
		b = a % b;
		a = temp;
	}
	return a;
}

function lcm(a, b) {
	if (a === b) {
		return a;
	} else {
		return (a * b) / gcd(a, b);
	}
}

function choices(n, m) {
	if (m === 0) {
		return [[]];
	} else if (m > n) {
		return [];
	} else if (m === n) {
		const choice = new Array(n);
		for (let i = 0; i < n; i++) {
			choice[i] = i;
		}
		return [choice];
	}
	const choicesIn = choices(n - 1, m - 1);
	const allChoices = choices(n - 1, m);
	for (let choice of choicesIn) {
		choice.push(n - 1);
		allChoices.push(choice);
	}
	return allChoices;
}

// Geometry calculations

/**
 * Line 1: An edge
 * Line 2: Perpendicular to the edge and running through another specified point
 * a1 * x + b1 * y + c1 = 0
 * a2 * x + b2 * y + c2 = 0
 * @return {number[]} [x2 - x1, a1, b1, c1, a2, b2]
 */
function projectionOntoLine(x1, y1, x2, y2) {
	const deltaX = x2 - x1;
	const deltaY = y2 - y1;
	if (deltaX === 0) {
		// x = -c1
		// y = -c2
		return [deltaX, 1, 0, -x1, 0, 1];

	} else if (deltaY === 0) {
		// y = -c1
		// x = -c2
		return [deltaX, 0, 1, -y1, 1, 0];
	}

	const a1 = -deltaY / deltaX;
	return [deltaX, a1, 1, -a1 * x1 - y1, deltaX / deltaY, 1];
}

function projectOntoLine(pointX, pointY, a1, b1, c1, a2, b2) {
	const c2 = -a2 * pointX - b2 * pointY;
	const denominator = a1 * b2 - a2 * b1;
	const intersectX = (b1 * c2 - b2 * c1) / denominator;
	const intersectY = (c1 * a2 - c2 * a1) / denominator;
	return [intersectX, intersectY];
}

function constrainToLineSegment(pointX, pointY, x1, y1, x2, y2) {
	const candidateX = [pointX, x1, x2];
	const candidateY = [pointY, y1, y2];
	let candidate = 0;
	if (x1 < x2) {
		if (pointX > x2) {
			candidate = 2;
		} else if (pointX < x1) {
			candidate = 1;
		}
	} else if (x2 < x1) {
		if (pointX > x1) {
			candidate = 1;
		} else if (pointX < x2) {
			candidate = 2;
		}
	} else {
		if (y1 < y2) {
			if (pointY > y2) {
				candidate = 2;
			} else if (pointY < y1) {
				candidate = 1;
			}
		} else if (y2 < y1) {
			if (pointY > y1) {
				candidate = 1;
			} else if (pointY < y2) {
				candidate = 2;
			}
		}
	}
	return [candidateX[candidate], candidateY[candidate]];
}

/**Finds contiguous runs of points which all share the same x-coordinate and sorts any such
 * subsequence by the y-coordinates of the points it contains.
 */
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
			subPointsY.sort(Sort.numericAscending);
		} else {
			subPointsY.sort(Sort.numericDescending);
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

function triangleArea(x1, y1, x2, y2, x3, y3) {
	return 0.5 * ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));
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
		const doubleTriArea = 2 * triangleArea(x0, y0, x1, y1, x2, y2);
		centreX += doubleTriArea * (x0 + x1 + x2);
		centreY += doubleTriArea * (y0 + y1 + y2);
		doubleArea += doubleTriArea;
	}
	centreX /= 3 * doubleArea;
	centreY /= 3 * doubleArea;
	return [centreX, centreY];
}

const Geometry = {
	choices: choices,
	gcd: gcd,
	lcm: lcm,
	centroid: centroid,
	constrainToLineSegment: constrainToLineSegment,
	orderPolygonPoints: orderPolygonPoints,
	projectionOntoLine: projectionOntoLine,
	projectOntoLine: projectOntoLine,
	triangleArea: triangleArea,
}

export {Geometry, Sort};
