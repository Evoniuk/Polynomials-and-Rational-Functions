class Polynomial:

    def __init__(self, coeff):
        self.coeff = coeff[::-1]
        while self.coeff[-1] == 0 and len(self.coeff) > 1:
            del self.coeff[-1]

    def __repr__(self):
        if len(self.coeff) == 1:
            return str(self.coeff[0])
        return 'Polynomial(' + str(self.coeff[::-1]) + ')'

    def __str__(self):
        if len(self.coeff) == 1:
            return str(self.coeff[0])
        return str(self.coeff[::-1])

    def __len__(self):
        return len(self.simplify().coeff) - 1

    def __eq__(self, other):
        if isinstance(other, (int, float)):
            if len(self) == 0 and self.coeff[0] == other:
                return True
            return False

        if isinstance(other, RationalFunction):
            return other == self

        def simplify(coeff):
            while coeff[-1] == 0 and len(coeff) > 1:
                del coeff[-1]
            return coeff

        return simplify(self.coeff[:]) == simplify(other.coeff[:])

    def __lt__(self, other):
        if self == other:
            return False
        if len(self.coeff) == len(other.coeff):
            for coeff in range(1, len(self.coeff) + 1):
                if self.coeff[-coeff] < other.coeff[-coeff]:
                    return True
                if self.coeff[-coeff] > other.coeff[-coeff]:
                    return False
        else:
            return len(self.coeff) < len(other.coeff)

    def __le__(self, other):
        return self < other or self == other

    def __add__(self, other):
        if isinstance(other, (int, float)):
            return self + Polynomial([other])
        if type(other) is RationalFunction:
            return other + self
        if len(self.coeff) > len(other.coeff):
            big = Polynomial(self.coeff[::-1])
            small = Polynomial(other.coeff[::-1])
        else:
            big = Polynomial(other.coeff[::-1])
            small = Polynomial(self.coeff[::-1])
        for a in range(len(small.coeff)):
            big.coeff[a] += small.coeff[a]
        return big

    def __neg__(self):
        return Polynomial([-coeff for coeff in self.coeff][::-1])

    def __sub__(self, other):
        return self + -other

    def __mul__(self, other):
        if isinstance(other, (int, float)):
            return self * Polynomial([other])
        if type(other) is RationalFunction:
            return other * self
        result = []
        polyResult = Polynomial([0])
        for coeff1 in range(len(self.coeff)):
            a = []
            for coeff2 in range(len(other.coeff)):
                a.append(self.coeff[coeff1] * other.coeff[coeff2])
            a = a[::-1]
            for b in range(coeff1):
                a.append(0)
            result.append(Polynomial(a))
        for polynomials in result:
            polyResult += polynomials
        return polyResult

    def __floordiv__(self, other):
        if type(other) is RationalFunction:
            if other.denom == Polynomial([1]):
                return self // other.numer
            else:
                raise ValueError("Polynomials are incompatible with rational functions.")
        if len(other) > len(self):
            return Polynomial([0])

        resultcoeff = self.coeff[-1] / other.coeff[-1]
        result = Polynomial([resultcoeff] + [0 for a in range(len(self) - len(other))])
        lower = self - result * other
        if lower == Polynomial([0]):
            return result
        return result + lower // other

    def __truediv__(self, other):
        if isinstance(other, (int, float)):
            return self * 1 / other
        if type(other) is RationalFunction:
            return self * other.reciprocal()
        return RationalFunction(self, other)

    def __mod__(self, other):
        if type(other) is RationalFunction:
            if other.denom == Polynomial([1]):
                return self % other.numer
            else:
                raise ValueError('Polynomials are incompatible with rational functions.')
        return (self - (self // other) * other).simplify()

    def simplify(self):
        while self.coeff[-1] == 0 and len(self.coeff) > 1:
            del self.coeff[-1]
        return self

    def evaluate(self, x):
        result = 0
        degree = 0
        for coeff in self.coeff:
            result += coeff * x**degree
            degree += 1
        return result

    def deriv_p(self):
        result = Polynomial([coeff for coeff in self.coeff[::-1]])
        for coeff in range(len(self.coeff)):
            result.coeff[coeff] *= coeff
        del result.coeff[0]
        return result

    def integrate(self):
        result = Polynomial([coeff for coeff in self.coeff[::-1]])
        result.coeff.insert(0, 0)
        for coeff in range(1, len(result.coeff)):
            result.coeff[coeff] /= coeff
        return result

    def slope(self, x):
        return self.deriv_p().evaluate(x)

    def concavity(self, x):
        if self.deriv_p().deriv_p().evaluate(x) > 0:
            return 'Concave Up'
        if self.deriv_p().deriv_p().evaluate(x) < 0:
            return 'Concave Down'
        return 'Inflection Point'

    def is_zero(self, x):
        return self.evaluate(x) == 0

    def factors(self):
        if len(self) == 1:
            return self
        if len(self) != 2:
            return None
        a = self.coeff[2]
        b = self.coeff[1]
        c = self.coeff[0]
        solution1 = (-b + (b ** 2 - 4 * a * c) ** (1 / 2)) / (2 * a)
        solution2 = (-b - (b ** 2 - 4 * a * c) ** (1 / 2)) / (2 * a)
        return (Polynomial([1, -solution1]), Polynomial([1, -solution2]))

    def int_factors(self):
        intfactors = set()
        if self.coeff[0] == 0:
            intfactors.add(0)
            tryagain = Polynomial(self.coeff[-1:0:-1])
            for x in tryagain.int_factors():
                intfactors.add(x)
        for x in range(abs(int(self.coeff[0])) + 1):
            denom = Polynomial([1, -x])
            if type((self / denom).simplify()) is Polynomial:
                intfactors.add(x)
            denom = Polynomial([1, x])
            if type((self / denom).simplify()) is Polynomial:
                intfactors.add(-x)
        return intfactors


class RationalFunction:

    def __init__(self, numer, denom):
        self.numer = numer
        self.denom = denom

    def __repr__(self):
        return 'RationalFunction(' + str(self.numer) + ', ' + str(self.denom) + ')'

    def __str__(self):
        frac = ''
        for a in range(max(len(str(self.numer)), len(str(self.denom)))):
            frac += '-'
        return str(self.numer) + '\n' + frac + '\n' + str(self.denom)

    def __eq__(self, other):
        if isinstance(other, (int, float)):
            return self == RationalFunction(Polynomial([other]), Polynomial([1]))
        if type(other) is Polynomial:
            return self == RationalFunction(other, Polynomial([1]))
        return self.numer == other.numer and self.denom == other.denom

    def __add__(self, other):
        if type(other) is Polynomial:
            return self + RationalFunction(other, Polynomial([1]))
        if isinstance(other, (int, float)):
            return self + Polynomial([other])
        newNum = self.numer * other.denom + other.numer * self.denom
        newDenom = self.denom * other.denom
        return RationalFunction(newNum, newDenom)

    def __neg__(self):
        return RationalFunction(-self.numer, self.denom)

    def __sub__(self, other):
        return self + -other

    def __mul__(self, other):
        if isinstance(other, (int, float)):
            newNum = [other * coeff for coeff in self.numer.coeff[::-1]]
            newNum = Polynomial(newNum)
            return RationalFunction(newNum, self.denom)
        if type(other) is Polynomial:
            return self * RationalFunction(other, Polynomial([1]))
        newNum = self.numer * other.numer
        newDenom = self.denom * other.denom
        return RationalFunction(newNum, newDenom)

    def __truediv__(self, other):
        if isinstance(other, (int, float)):
            return self * (1 / other)
        if type(other) is Polynomial:
            return self / RationalFunction(other, Polynomial([1]))
        return self * other.reciprocal()

    def simplify(self):
        quotient = self.numer // self.denom
        remainder = self.numer % self.denom
        if remainder == Polynomial([0]):
            return quotient
        return self

    def reciprocal(self):
        return RationalFunction(self.denom, self.numer)

    def evaluate(self, x):
        num = 0
        denom = 0
        for coeff in range(len(self.numer.coeff)):
            num += self.numer.coeff[coeff] * x**coeff
        for coeff in range(len(self.denom.coeff)):
            denom += self.denom.coeff[coeff] * x**coeff
        return num / denom

    def deriv_r(self):
        result = RationalFunction(self.numer, self.denom)
        result.numer = result.denom * result.numer.deriv_p() - result.numer * result.denom.deriv_p()
        result.denom *= result.denom
        return result

    def slope(self, x):
        return self.deriv_r().evaluate(x)

    def concavity(self, x):
        if self.deriv_r().deriv_r().evaluate(x) > 0:
            return 'Concave Up'
        if self.deriv_r().deriv_r().evaluate(x) < 0:
            return 'Concave Down'
        return 'Inflection Point'

    def is_zero(self, x):
        return self.evaluate(x) == 0


# poly1 = Polynomial([1, 7, 12])
# poly2 = Polynomial([1, 3])
# ratio = RationalFunction(poly1, poly2)

# print(poly1)
# print(poly2)
# print(ratio)
# print(ratio.simplify())
