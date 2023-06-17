class Shape {

	constructor(pointsX, pointsY) {
		const numPoints = pointsX.length;
		let sumX = 0, sumY = 0;
		for (let i = 0; i < numPoints; i++) {
			sumX += pointsX[i];
			sumY += pointsY[i];
		}
		const centreX = sumX / numPoints;
		const centreY = sumY / numPoints;

		const deltaX = new Array(numPoints);
		const deltaY = new Array(numPoints);
		for (let i = 0; i < numPoints; i++) {
			deltaX[i] = pointsX[i] - centreX;
			deltaY[i] = pointsY[i] - centreY;
		}

		const angles = new Array(numPoints);
		const indexes = new Array(numPoints);
		for (let i = 0; i < numPoints; i++) {
			angles[i] = Math.atan2(deltaY[i], deltaX[i]);
			indexes[i] = i;
		}

		indexes.sort((i, j) => angles[i] - angles[j]);

		const sortedX = new Array(numPoints);
		const sortedY = new Array(numPoints);
		const sortedDeltaX = new Array(numPoints);
		const sortedDeltaY = new Array(numPoints);;
		const sortedAngles = new Array(numPoints);
		const radii = new Array(numPoints);
		let sumR = 0;
		for (let i = 0; i < numPoints; i++) {
			const index = indexes[i];
			sortedX[i] = pointsX[index];
			sortedY[i] = pointsY[index];
			const dx = deltaX[index];
			const dy = deltaY[index];
			sortedDeltaX[i] = dx;
			sortedDeltaY[i] = dy;
			sortedAngles[i] = angles[index];
			const radius = Math.hypot(dx, dy);
			radii[i] = radius;
			sumR += radius;
		}

		this.numPoints = numPoints;
		this.pointsX = sortedX;
		this.pointsY = sortedY;
		this.centreX = centreX;
		this.centreY = centreY;
		this.deltaX = sortedDeltaX;
		this.deltaY = sortedDeltaY;
		this.angles = sortedAngles;
		this.radii = radii;
		this.size = sumR / numPoints;

		this.resizedX = sortedDeltaX;
		this.resizedY = sortedDeltaY;
		this.rotatedX = sortedDeltaX;
		this.rotatedY = sortedDeltaY;
		this.rotation = 0;
		// Align Vertex 0 of this shape with this numbered vertex in the target shape.
		this.vertexRotations = 0;
		// The offsets that remain after translation, resizing and rotating.
		this.offsetsX = undefined;
		this.offsetsY = undefined;
	}

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

	rotate(shape2) {
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
			if (Math.abs(rotation) > 0.5 * Math.PI) {
				continue;
			}
			const cos = Math.cos(rotation);
			const sin = Math.sin(rotation);
			const rotatedX = new Array(numPoints);
			const rotatedY = new Array(numPoints);
			let distanceSquared = 0;
			for (let j = 1; j < numPoints; j++) {
				const x = this.resizedX[j];
				const y = this.resizedY[j];
				const transformedX = x * cos - y * sin;
				const transformedY = x * sin + y * cos;
				const vertexIndex2 = (i + j) % numPoints2;
				const targetX = shape2.resizedX[vertexIndex2];
				const targetY = shape2.resizedY[vertexIndex2];
				const distanceX = targetX - transformedX;
				const distanceY = targetY - transformedY;
				distanceSquared += distanceX * distanceX + distanceY * distanceY;
				rotatedX[j] = transformedX;
				rotatedY[j] = transformedY;
			}
			if (distanceSquared < minError) {
				minError = distanceSquared;
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
		pointsX[i] = Math.random() * canvas.width;
		pointsY[i] = Math.random() * canvas.height;
	}
	return new Shape(pointsX, pointsY);
}

const canvas = document.getElementById('canvas');
const context = canvas.getContext('2d');
canvas.width = window.innerWidth;
canvas.height = window.innerHeight;

function drawPolygon(pointsX, pointsY) {
	const numPoints = pointsX.length;
	context.beginPath();
	context.moveTo(pointsX[0], pointsY[0]);
	for (let i = 1; i < numPoints; i++) {
		const x = pointsX[i];
		const y = pointsY[i];
		context.lineTo(x, y);
	}
	context.fill();
}

function findInterpolationStep(polygon1, polygon2, speed = 1) {
	const numPoints = polygon2.numPoints;
	const inverseScale = polygon1.size / polygon2.size;

	const translateX = polygon2.centreX - polygon1.centreX;
	const translateY = polygon2.centreY - polygon1.centreY;
	let maxTranslation = 0, maxArcLength = 0;
	for (let i = 0; i < numPoints; i++) {
		maxTranslation = Math.max(
			maxTranslation,
			Math.abs(translateX + polygon2.offsetsX[i]),
			Math.abs(translateY + polygon2.offsetsY[i])
		);

		const radius = polygon2.radii[i];
		const scaledRadius = radius * inverseScale;
		// Arc length of an Archimedean spiral.
		const radialGrowth = (scaledRadius - radius) / Math.abs(polygon2.rotation);
		let squareRoot = Math.sqrt(scaledRadius * scaledRadius + radialGrowth * radialGrowth);
		const firstTerm =  0.5 * squareRoot * scaledRadius / radialGrowth;
		const secondTerm = 0.5 * radialGrowth * Math.log(squareRoot + scaledRadius);
		squareRoot = Math.sqrt(radius * radius + radialGrowth * radialGrowth);
		const thirdTerm = 0.5 * squareRoot * radius / radialGrowth;
		const fourthTerm = 0.5 * radialGrowth * Math.log(squareRoot + radius);
		const arcLength = (firstTerm + secondTerm) - (thirdTerm + fourthTerm);
		maxArcLength = Math.max(maxArcLength, arcLength);
	}
	return speed / Math.max(maxTranslation, maxArcLength);
}

function draw(interpolation) {
	context.clearRect(0, 0, canvas.width, canvas.height);

	context.globalAlpha = 0.4;
	context.fillStyle = 'red';
	drawPolygon(polygon1.pointsX, polygon1.pointsY);
	context.fillStyle = 'blue';
	drawPolygon(polygon2.pointsX, polygon2.pointsY);

	let translateX, translateY, scale, rotation;
	translateX = polygon2.centreX * interpolation + polygon1.centreX * (1 - interpolation);
	translateY = polygon2.centreY * interpolation + polygon1.centreY * (1 - interpolation);
	scale = polygon2.size / polygon1.size * interpolation + 1 - interpolation;
	rotation = -polygon2.rotation * interpolation;

	const numPoints = polygon2.numPoints;
	const pointsX = new Array(numPoints);
	const pointsY = new Array(numPoints);
	for (let i = 0; i < numPoints; i++) {
		pointsX[i] = polygon2.rotatedX[i] - (1 - interpolation) * polygon2.offsetsX[i];
		pointsY[i] = polygon2.rotatedY[i] - (1 - interpolation) * polygon2.offsetsY[i];
	}

	context.fillStyle = 'lime';
	context.save();
	context.translate(translateX, translateY);
	context.scale(scale, scale);
	context.rotate(rotation);
	drawPolygon(pointsX, pointsY);
	context.restore();
}

let interpolationStep, animationStart;

function animate(time) {
	let interpolation;
	if (animationStart === undefined) {
		animationStart = time;
		interpolation = 0;
	} else {
		interpolation = Math.min((time - animationStart) * interpolationStep, 1);
	}
	draw(interpolation);
	if (interpolation < 1) {
		requestAnimationFrame(animate);
	}
}

const polygon1 = randomPolygon(5);
const polygon2 = randomPolygon(5);
polygon2.resize(polygon1.size);
polygon2.rotate(polygon1);
interpolationStep = findInterpolationStep(polygon1, polygon2, 0.25);
requestAnimationFrame(animate);
