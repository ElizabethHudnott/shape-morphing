import {Geometry} from './math.js';

function tracePolygonPath(context, pointsX, pointsY) {
	const numPoints = pointsX.length;
	context.beginPath();
	context.moveTo(pointsX[0], pointsY[0]);
	for (let i = 1; i < numPoints; i++) {
		const x = pointsX[i];
		const y = pointsY[i];
		context.lineTo(x, y);
	}
}

function tracePolylinePath(context, pointsX, pointsY, closed = true, start = 0, end = 1) {
	const segmentOffsets = Geometry.getPerimeterOffsets(pointsX, pointsY, closed);
	const [startVertex, startProportion, startX, startY] = Geometry.getPerimeterOffset(
		start, pointsX, pointsY, segmentOffsets, closed
	);
	const [endVertex, endProportion, endX, endY] = Geometry.getPerimeterOffset(
		end, pointsX, pointsY, segmentOffsets, closed
	);
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

function drawInterpolatedShape(context, morph, interpolation) {
	tracePolygonPath(context, morph.pointsX, morph.pointsY);
	context.fill();

	const strokeMorph = morph.strokeMorph;
	if (strokeMorph?.width > 0) {
		tracePolylinePath(
			context, morph.pointsX, morph.pointsY, strokeMorph.closed, strokeMorph.start,
			strokeMorph.end
		);
		context.globalAlpha = 1;
		context.shadowColor = 'transparent';
		context.stroke();
	}
}

function drawVertexMap(context, polygon1, polygon2) {
	const numPoints = polygon1.numPoints;
	context.clearRect(0, 0, canvas.width, canvas.height);
	context.save();
	context.translate(0.5 * canvas.width,  0.5 * canvas.height);

	// Draw first polygon in red.
	context.beginPath();
	context.moveTo(polygon1.deltaX[0], polygon1.deltaY[0]);
	for (let i = 1; i < numPoints; i++) {
		context.lineTo(polygon1.deltaX[i], polygon1.deltaY[i]);
	}
	context.closePath();
	context.fillStyle = 'rgba(255, 0, 0, 70%)';
	context.fill();

	// Draw second polygon in green.
	context.beginPath();
	context.moveTo(polygon2.rotatedX[0], polygon2.rotatedY[0]);
	for (let i = 1; i < numPoints; i++) {
		context.lineTo(polygon2.rotatedX[i], polygon2.rotatedY[i]);
	}
	context.closePath();
	context.fillStyle = 'rgba(0, 255, 0, 70%)';
	context.fill();

	// Draw lines between equal numbered vertices.
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


export default {
	tracePolygonPath: tracePolygonPath,
	tracePolylinePath: tracePolylinePath,
	drawInterpolatedShape: drawInterpolatedShape,
	drawVertexMap: drawVertexMap,
};
