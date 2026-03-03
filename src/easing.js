// Easing functions.

// Circular shape.

function circleIn(m = Infinity) {
	if (m === Infinity) {
		return t => 1 - Math.sqrt(1 - t * t);
	}
	const x = 0.5 * (Math.sqrt(4 * m * m + 1) - 1) / m;
	return t => (1 - Math.sqrt(1 - x * x * t * t)) / x;
}

function circleOut(m = Infinity) {
	const f = circleIn(m);
	return t => f(1) - f(1 - t)
};

function circleInOut(m = Infinity, x = 0.5) {
	if (m === Infinity) {
		return twoPart(circleIn(), circleOut(), x);
	}
	const xm = 0.5 * (Math.sqrt(4 * m * m + 1) - 1) / m;
	const bigCircle = t => (1 - Math.sqrt(1 - xm * xm * t * t)) / xm;
	const bigCircle1 = bigCircle(1);
	const out = t => bigCircle1 - bigCircle(1 - t);
	const unscaled =  t => t <= x ? bigCircle(t / x) : out((t - x) / x) + bigCircle1;
	const max = unscaled(1);
	return t => unscaled(t) / max;
}

// Exponential

function expIn(k) {
	return t => -Math.exp(-k) * Math.expm1(k * t) / Math.expm1(-k);
}

function expOut(k) {
	return t => Math.expm1(-k * t) / Math.expm1(-k);
}

const expInOut = k => twoPart(expIn(k), expOut(k));

// Linear

const linear = t => t;

// Quadratic

const quadIn =  t => t * t;
const quadOut = t => t * (2 - t);
const quadInOut = twoPart(quadIn, quadOut);

// Quarter and half sine shapes

const sineIn  = t => 1 - Math.cos(0.5 * Math.PI * t);
const sineOut = t => Math.sin(0.5 * Math.PI * t);
const sineInOut = t => 0.5 - 0.5 * Math.cos(Math.PI * t);

// Step functions

function stepsJumpStart(n) {
	return t => Math.ceil(t * n) / n;
}

function stepsJumpEnd(n) {
	return t => Math.trunc(t * n) / n;
}

// Ways of combining multiple easing functions.

function twoPart(f, g, x = 0.5, y = 0.5) {
	return t => t <= x ? y * f(t / x) : (1 - y) * g((t - x) / (1 - x)) + y;
}

function threePart(f, g, h, x1, y1, x2, y2) {
	function ease(t) {
		if (t <= x1) {
			return y1 * f(t / x1);
		} else if (t <= x2) {
			return (y2 - y1) * g((t - x1) / (x2 - x1)) + y1;
		} else {
			return (1 - y2) * h((t - x2) / (1 - x2)) + y2;
		}
	}
	return ease;
}

export default {
	circleIn: circleIn,
	circleOut: circleOut,
	circleInOut: circleInOut,
	expIn,
	expOut,
	expInOut,
	linear,
	quadIn,
	quadOut,
	quadInOut,
	sineIn,
	sineOut,
	sineInOut,
	stepsJumpStart,
	stepsJumpEnd,
	twoPart,
	threePart,
};
