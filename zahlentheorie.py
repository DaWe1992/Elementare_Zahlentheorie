# -*- coding: utf-8 -*-
"""
Spyder Editor

Algorithms Number Theory
"""
import math

from itertools import combinations
from functools import reduce
from operator import mul


def gcd(a, b, p=True):
    """
    Finds greates common divisor of a and b.
    
    :param a:           first number
    :param b:           second number
    :param p:           flag indicating whether to print steps
    :return:            greatest common divisor
    """
    a = abs(a)
    b = abs(b)
    
    while b != 0:
        q = a // b
        r = a % b
        if p:
            print("{:03d} = {:03d} * {:03d} + {:03d}".format(a, q, b, r))
        a = b
        b = r
        
    return a


def gcdExtended(a, b):
    """
    Extended euclidean algorithm.
    
    :param a:           first number
    :param b:           second number
    :return:            greatest common divisor,
                        s, t such that gcd(a, b) = a * s + b * t
    """
    a = abs(a)
    b = abs(b)
    
    if a == 0 :  
        return b, 0, 1
             
    gcd, x1, y1 = gcdExtended(b % a, a) 
    
    x = y1 - (b // a) * x1 
    y = x1 
     
    return gcd, x, y


def lcm(a, b):
    """
    Finds least common multiple of a and b
    
    :param a:           first number
    :param b:           second number
    :return:            least common multiple
    """
    return (abs(a) * abs(b) / gcd(a, b, False))


def diophant(a, b, c):
    """
    Solves linear diophant equation aX + bY = c.
    
    :param a:           coefficient a
    :param b:           coefficient b
    :param c:           coefficient c
    """
    gcd, s, t = gcdExtended(a, b)
    print("ggT(a, b) = {}".format(gcd))
    if c % gcd == 0:
        print("ggT(a, b) = {} | {}".format(gcd, c))
        print("Es existieren unendlich viele Lösungen")
        m = c / gcd
        print("m = c / ggT(a, b) = {} / {} = {}".format(c, gcd, m))
        print("Es ist ggT({0}, {1}) = {0} * {2} + {1} * {3} = {4}".format(a, b, s, t, gcd))
        x0 = m * s
        y0 = m * t
        print("x0 = {} * {} = {}".format(m, s, x0))
        print("y0 = {} * {} = {}".format(m, t, y0))
        a_ = a / gcd
        b_ = b / gcd
        print("a' = {} / {} = {}".format(a, gcd, a_))
        print("b' = {} / {} = {}".format(b, gcd, b_))
        
        print("L = {{ ({0} + {1} * k, {2} - {3} * k) | k in Z }}".format(x0, b_, y0, a_))
    else:
        print("Es existiert keine Lösung")
        
        
def congr(a, b, m, p=True):
    """
    Solves linear congruence.
    
    :param a:           number 1
    :param b:           number 2
    :param c:           modulus
    :param p:           flag indicating whether to print steps
    :return:            solutions
    """
    sols = []
    d = gcd(a, m, False)
    if p: print("ggT(a, m) = {}".format(d))
    
    if b % d == 0:
        if p: print("ggT(a, m) = {} | {}".format(d, b))
        if p: print("Es existieren {} Lösungen".format(d))
        _, s, t = gcdExtended(a / d, m / d)
        if p: print("1 = {} * {} + {} * {}".format(a / d, s, m / d, t))
        x0 = s * (b / d) % (m / d)
        if p: print("x0 = s * b / d mod m / d = {0} * {1} / {2} mod {3} / {2} = {4}".format(s, b, d, m, x0))
        sols.append(x0)
        
        for i in range(1, d):
            sol = x0 + i * m / d
            if p: print("x{0} = {1} + {0} * {2} / {3} = {4}".format(i, x0, m, d, sol))
            sols.append(sol)
            
        if p: print(sols)
        return sols
    else:
        print("Es existiert keine Lösung")
        
        
def chineseRem(a_list, m_list):
    """
    Chinese remainder.
    
    :param a_list:      list remainder
    :param m_list:      list modulus
    """
    for pair in list(combinations(m_list, 2)):
        if gcd(pair[0], pair[1], False) > 1:
            print("Module nicht paarweise teilerfremd")
            return
    
    # calculate M
    M = reduce(mul, m_list, 1)
    print("M = {}".format(M))
    
    res = 0
    
    for i in range(len(m_list)):
        # calculate Mi
        Mi = int(M / m_list[i])
        print("M{} = {}".format(i + 1, Mi))
        # calculate xi
        x = congr(Mi, 1, m_list[i], False)[0]
        print("{}X = 1 (mod {}) => x{} = {}".format(Mi, m_list[i], i + 1, x))
        
        res += a_list[i] * Mi * x
        
    print("Lösung: {}".format(res % M))
    
    
def sumOfSquares(p):
    """
    Calculates a and b such that p = a^2 + b^2 (if possible)
    
    :param p:           prime number
    :return:            [a, b]
    """
    def KtoK1(a, b, k):
        print("\nSchritt 2:\n==========")
        # find r
        r = congr(1, a, k, False)[0]
        print("Wir teilen a durch k mit Rest:")
        print("Es ist r~ = {}".format(r))
        if r > k / 2:
            print("r~ > {0} / 2 -> r = {1} - {0}".format(k, r))
            r = r - k
        else:
            print("-{0}/2 < {1} <= {0}/2 -> Setze r = {1}".format(k, r))
            
        # find s
        s = congr(1, b, k, False)[0]
        print("Wir teilen b durch k mit Rest:")
        print("Es ist s~ = {}".format(s))
        if s > k / 2:
            print("s~ > {0} / 2 -> s = {1} - {0}".format(k, s))
            s = s - k
        else:
            print("-{0}/2 < {1} <= {0}/2 -> Setze s = {1}".format(k, s))
            
        print("\nSchritt 3:\n==========")
        print("r und s sind nicht beide 0")
        
        print("\nSchritt 4:\n==========")
        k1 = int((r**2 + s**2) // k)
        print("k1 = (r^2 + s^2) / k = ({0}^2 + {1}^2) / {2} = {3}".format(r, s, k, k1))
    
        print("\nSchritt 5:\n==========")
        print("""
              (ra + sb)^2 + (rb - sa)^2 = ({0} * {1} + {2} * {3})^2 + ({0} * {3} - {2} * {1})^2
                                        = ({4})^2 + ({5})^2
                                        = {6} * {7}^2 * {8} = {9}
                                        = k1 * k^2 * p
        """.format(r, a, s, b, r * a + s * b, r * b - s * a, k1, k, p, k1 * k**2 * p))
        
        print("\nSchritt 6:\n==========")
        a_ = int((r * a + s * b) // k)
        print("a = (r * a + s * b) / k = ({0} * {1} + {2} * {3}) / {4} = {5}".format(r, a, s, b, k, a_))
        b_ = int((r * b - s * a) // k)
        print("b = (r * b - s * a) / k = ({0} * {1} - {2} * {3}) / {4} = {5}".format(r, b, s, a, k, b_))
        return abs(a_), abs(b_), k1
    
    
    print("Schritt 1:\n==========")
    if p % 4 != 1:
        print("{} kann nicht als Summe von zwei Quadraten geschrieben werden.".format(p))
        return []
    print("{} ist kongruent 1 modulo 4".format(p))
    a = math.factorial((p - 1) / 2) % p
    print("a = (({0} - 1) / 2)! mod {0} = {1}".format(p, a))
    b = 1
    print("b = 1")
    k = (a**2 + b**2) // p
    print("k = ({0}^2 + {1}^2) / {2} = {3}".format(a, b, p, k))
    
    while k > 1:
        a, b, k = KtoK1(a, b, k)
        
    print("\nErgebnis:")
    print("{0} = {1}^2 + {2}^2".format(p, a, b))
    return [a, b]


def sumOf4Squares(p):
    """
    Calculates sum of four squares.
    
    :param p:       prime number
    :return:        [a, b, c, d]
    """
    
    def lemma4quads(p):
        """
        Lemma 6.3.4
        """
        print("\nLemma 6.3.4:")
        k = ((p - 1) // 2)
        print("({0} - 1) / 2 = {1}".format(p, k))
        S1 = []
        S2 = []
        
        for i in range(k + 1):
            S1.append(i**2 % p)
            S2.append((-1 - i**2) % p)
            
        print("Wir betrachten die Menge S1 = {{0^2 mod {0}, ..., {1}^2 mod {0}}}".format(p, k))
        print("S1 =", S1)
        print("Wir betrachten die Menge S2 = {{-1-0^2 mod {0}, ..., -1-{1}^2 mod {0}}}".format(p, k))
        print("S2 =", S2)
        
        for i in range(k + 1):
            for j in range(k + 1):
                if S1[i] == S2[j]:
                    print("Es ist 1 + {0}^2 + {1}^2 kongruent 0 (mod {2})".format(i, j, p))
                    return i, j
                
                
    def mUngerade(a, b, c, d, p):
        """
        Lemma 6.3.7
        """
        print("\nLemma 6.3.7:")
        print("Umsortieren der Summanden, sodass x mod 2 = y mod 2 und z mod 2 = w mod 2:")
        if a % 2 == b % 2:
            x = a; y = b; z = c; w = d
        elif a % 2 == c % 2:
            x = a; y = c; z = b; w = d
        else:
            x = a; y = d; z = c; w = b
        print("x = {0}, y = {1}, z = {2}, w = {3}\n".format(x, y, z, w))
            
        k = (a**2 + b**2 + c**2 + d**2) // p
        
        while k % 2 == 0:
            print("k = ({0}^2 + {1}^2 + {2}^2 + {3}^2) / {4} = {5}".format(x, y, z, w, p, k))
            if k % 2 == 0: print("k ist gerade")
            else: print("k ist ungerade")
        
            x1 = (x - y) // 2
            y1 = (x + y) // 2
            z1 = (z - w) // 2
            w1 = (z + w) // 2
            print("""
                x1 = ({0} - {1}) / 2 = {4}
                y1 = ({0} + {1}) / 2 = {5}
                z1 = ({2} - {3}) / 2 = {6}
                w1 = ({2} + {3}) / 2 = {7}
            """.format(x, y, z, w, x1, y1, z1, w1))
            print("({0} / 2) * {1} = ({2})^2 + ({3})^2 + ({4})^2 + ({5})^2\n"
                  .format(k, p, x1, y1, z1, w1))
            
            if x1 % 2 == y1 % 2:
                x = x1; y = y1; z = z1; w = w1
            else:
                x = x1; y = w1; z = z1; w = y1
                
            k = (x**2 + y**2 + z**2 + w**2) // p
          
        print("k = (({0})^2 + ({1})^2 + ({2})^2 + ({3})^2) / {4} = {5}".format(x, y, z, w, p, k))
        print("k ist ungerade")
        print("Ergebnis:", abs(x), abs(y), abs(z), abs(w))
        return abs(x), abs(y), abs(z), abs(w), k
    
    
    def k2k1(a, b, c, d, p):
        """
        Lemma 6.3.10
        """
        print("\nLemma 6.3.10:")
        k = (a**2 + b**2 + c**2 + d**2) // p
        print("k = {} > 1".format(k))
        r = mods(a, k)
        s = mods(b, k)
        t = mods(c, k)
        u = mods(d, k)
        print("""
            r = mods(a, k) = mods({0}, {4}) = {5}
            s = mods(b, k) = mods({1}, {4}) = {6}
            t = mods(c, k) = mods({2}, {4}) = {7}
            u = mods(d, k) = mods({3}, {4}) = {8}
        """.format(a, b, c, d, k, r, s, t, u))
        
        a1 = abs((a * r + b * s + c * t + d * u) // k)
        b1 = abs((a * s - b * r + c * u - d * t) // k)
        c1 = abs((a * t - b * u - c * r + d * s) // k)
        d1 = abs((a * u + b * t - c * s - d * r) // k)
        print("""
            a1 = |(a * r + b * s + c * t + d * u) / k| = |{0} * {1} + {2} * {3} + {4} * {5} + {6} * {7}) / {8}| = {9}
            b1 = |(a * s - b * r + c * u - d * t) / k| = |{0} * {3} + {2} * {1} + {4} * {7} + {6} * {5}) / {8}| = {10}
            c1 = |(a * t - b * u - c * r + d * s) / k| = |{0} * {5} + {2} * {7} + {4} * {1} + {6} * {3}) / {8}| = {11}
            d1 = |(a * u + b * t - c * s - d * r) / k| = |{0} * {7} + {2} * {5} + {4} * {3} + {6} * {1}) / {8}| = {12}
        """.format(a, r, b, s, c, t, d, u, k, a1, b1, c1, d1))
        
        return mUngerade(a1, b1, c1, d1, p)
        
    
    print("Schreibe {} als Summe von vier Quadraten:".format(p))
    if p == 2: return 1, 1, 0, 0
    u, v = lemma4quads(p)
    print("u = {0}, v = {1}".format(u, v))
    print("Summanden: {0} {1} 1 0".format(u, v))
    a, b, c, d, k = mUngerade(u, v, 1, 0, p)

    while k > 1:
        a, b, c, d, k = k2k1(a, b, c, d, p)
    
    print("\nErgebnis:")
    print("{0} = {1}^2 + {2}^2 + {3}^2 + {4}^2".format(p, a, b, c, d))
    return [a, b, c, d]

    
def mods(a, b):
    """
    Calculates symmetric mod.
    
    :param a:       dividend
    :param b:       divisor
    :return:        symmetric mod
    """
    m = a % b
    if m > b / 2: m -= b
    return m


def pellEq(x1, y1, d, n, verbose=True):
    """
    Finds solutions for Pell's equation.
    
    :param x0:      initial solution (part 1)
    :param y0:      initial solution (part 2)
    :param d:       d parameter in pell's equation
    :param n:       number of solutions to calculate
    """
    print("1. Lösung:", (x1, y1), "\n")
    x_ = x1
    y_ = y1
    
    for i in range(n):
        x2 = x1 * x_ + d * y1 * y_
        y2 = x1 * y_ + x_ * y1
        if verbose: print("a{0} = {1} * {2} + {3} * {4} * {5} = {6}".format(i + 2, x1, x_, d, y1, y_, x2))
        if verbose: print("b{0} = {1} * {2} + {3} * {4} = {5}".format(i + 2, x1, y_, x_, y1, y2))
        print("{}. Lösung:".format(i + 2), (x2, y2), "\n")
        x_ = x2
        y_ = y2
        

if __name__ == "__main__":
    # print(gcd(20, 12))    
    # print(gcdExtended(20, 12))
    # print(lcm(9, 4))
    # diophant(2, 6, 20)
    # congr(246, 18, 303)
    # chineseRem([3, 2, 6], [8, 5, 9])
    # sumOfSquares(21) # 1277
    print(sumOf4Squares(23))
    # pellEq(3, 2, 2, 5, verbose=False)
    