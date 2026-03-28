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


export default {
	tracePolygonPath: tracePolygonPath,
	tracePolylinePath: tracePolylinePath,
	drawInterpolatedShape: drawInterpolatedShape,
	drawVertexMap: drawVertexMap,
};
