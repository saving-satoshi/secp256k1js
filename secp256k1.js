'use strict';

// Most of this code is based on:
// https://github.com/bitcoin/bitcoin/blob/1830dd8820fb90bac9aea32000e47d7eb1a99e1b/test/functional/test_framework/secp256k1.py
// Copyright (c) 2022-2023 The Bitcoin Core developers
// Distributed under the MIT software license, see the accompanying
// file COPYING or http://www.opensource.org/licenses/mit-license.php.

const {inspect} = require('util');
const FIELD_SIZE = 2n**256n - 2n**32n - 977n;
const ORDER = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141n;
const ORDER_HALF = ORDER / 2n;

class FE {
  constructor(val) {
    if (typeof val !== 'BigInt')
      val = BigInt(val);
    this.val = ((val % FIELD_SIZE) + FIELD_SIZE) % FIELD_SIZE;
  }

  add(n) {
    return new FE(this.val + (n.val || n));
  }

  mul(n) {
    return new FE(this.val * (n.val || n));
  }

  sub(n) {
    return new FE(this.val - (n.val || n));
  }

  div(n) {
    if (!n.val)
      n = new FE(n);

    return this.mul(n.inv());
  }

  pow(exponent) {
    if (exponent.val)
      exponent = exponent.val;

    let result = 1n;
    let base = this.val % FIELD_SIZE;
    while (exponent > 0n) {
      if (exponent % 2n === 1n)  // odd number
        result = (result * base) % FIELD_SIZE;
      exponent = exponent >> 1n; // divide by 2
      base = (base * base) % FIELD_SIZE;
    }
    return new FE(result);
  }

  // Modular multiplicative inverse using Extended Euclidean Algorithm
  inv() {
    let x0 = 0n;
    let x1 = 1n;
    let a = this.val;
    let m = FIELD_SIZE;

    while (a > 1n) {
      const q = a / m;
      let t = m;
      m = a % m;
      a = t;
      t = x0;
      x0 = x1 - q * x0;
      x1 = t;
    }

    if (x1 < 0n)
      x1 += FIELD_SIZE;

    return new FE(x1);
  }

  // Compute the square root of a field element if it exists (None otherwise).
  // Due to the fact that our modulus is of the form (p % 4) == 3, the Tonelli-Shanks
  // algorithm (https://en.wikipedia.org/wiki/Tonelli-Shanks_algorithm) is simply
  // raising the argument to the power (p + 1) / 4.
  // To see why: (p-1) % 2 = 0, so 2 divides the order of the multiplicative group,
  // and thus only half of the non-zero field elements are squares. An element a is
  // a (nonzero) square when Euler's criterion, a^((p-1)/2) = 1 (mod p), holds. We're
  // looking for x such that x^2 = a (mod p). Given a^((p-1)/2) = 1, that is equivalent
  // to x^2 = a^(1 + (p-1)/2) mod p. As (1 + (p-1)/2) is even, this is equivalent to
  // x = a^((1 + (p-1)/2)/2) mod p, or x = a^((p+1)/4) mod p.
  sqrt() {
    const s = this.pow((FIELD_SIZE + 1n) / 4n);
    if (s.pow(2n).val === this.val)
      return s
    return null;
  }

  equals(n) {
    return this.val === n.val;
  }

  neg() {
    return new FE(-this.val);
  }

  isEven() {
    return !(this.val & 1n)
  }

  hex() {
    return this.val.toString(16).padStart(64, '0');
  }

  getBytes() {
    return Buffer.from(this.hex(), 'hex');
  }

  [inspect.custom]() {
    return `<FE: ${this.hex()}>`;
  }
}

class GE {
  constructor(x, y) {
    this.inf = true;
    this.x = null;
    this.y = null;

    if (x && y) {
      if (!y.pow(2n).equals(x.pow(3n).add(7n)))
        throw new Error('Point is not on curve.')
      this.inf = false;
      this.x = x;
      this.y = y;
    }
  }

  add(a) {
    // Deal with infinity: a + infinity == infinity + a == a.
    if (this.inf)
      return a;
    if (a.inf)
      return this;

    let lam;
    if (this.x.equals(a.x)) {
      if (!this.y.equals(a.y)) {
        // A point added to its own negation is infinity.
        if (!this.y.add(a.y).equals(0n))
          throw new Error('One point is not on curve.')
        return new GE();
      } else {
        // For identical inputs, use the tangent (doubling formula).
        lam = this.x.pow(2n).mul(3n).div(this.y.mul(2n));
      }
    } else {
      // For distinct inputs, use the line through both points (adding formula).
      lam = this.y.sub(a.y).div(this.x.sub(a.x));
    }

    // Determine point opposite to the intersection of that line with the curve.
    const x = lam.pow(2n).sub(this.x.add(a.x));
    const y = lam.mul(this.x.sub(x)).sub(this.y);
    return new GE(x, y);
  }

  mul(a) {
    // Reduce scalar modulo order of the curve
    a = a % ORDER;
    // Start with the point at infinty
    let r = new GE();
    // Iterate over all bit positions from high to low
    for (let i = 255n; i >= 0n; i--) {
      // Double
      r = r.add(r);
      // Then add the points for which corresponding scalar bit is set
      if (a >> i & 1n)
        r = r.add(this)
    }

    return r;
  }

  static liftX(x) {
    let y = x.pow(3n).add(7n).sqrt();
    if (!y)
      return null;
    if (!y.isEven())
      y = y.neg();
    return new GE(x, y);
  }

  [inspect.custom]() {
    return `<GE: ${this.x.hex()}, ${this.y.hex()}>`;
  }
}

exports.FIELD_SIZE = FIELD_SIZE;
exports.FE = FE;
exports.GE = GE;
exports.G = GE.liftX(new FE(0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798n));
exports.ORDER = ORDER;
