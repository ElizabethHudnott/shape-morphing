import Ease from './easing.js';
import {Geometry} from './math.js';

class ConstantColour {

	static BLACK = new ConstantColour('black');

	constructor(colour) {
		this.colour = colour;
	}

	interpolate(morph, interpolation) {
		return this.colour;
	}
}

class Colour {

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

class VertexVertexGradient {

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

class EdgeParallelGradient {

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

class EdgePerpendicularGradient {

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

class VertexConicGradient {

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

class CentreConicGradient {

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


export default {
	Constant: ConstantColour,
	Colour: Colour,
	VertexVertexGradient: VertexVertexGradient,
	EdgeParallelGradient: EdgeParallelGradient,
	EdgePerpendicularGradient: EdgePerpendicularGradient,
	VertexConicGradient: VertexConicGradient,
	CentreConicGradient: CentreConicGradient,
};
