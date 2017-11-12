fun calcFunction(x: Double): Double {
    return 1 - Math.sin(x)
}

fun calcKernel(x: Double, s: Double): Double {
    return Math.cos(x - s)
}

fun gauss(aMtr: Matrix, bVect: Vector): Vector {
    val a = Matrix(aMtr)
    val b = Vector(bVect)
    val x = Vector(b.size)
    val xIndexes = Vector(b.size)
    var max: Double
    var maxk: Int
    var index: Int

    for (i in 0 until x.size) {
        xIndexes[i] = i.toDouble()
    }

    for (k in 0 until x.size) {
        max = a[k][k]
        maxk = k
        for (i in k + 1 until x.size) {
            if (Math.abs(max) < Math.abs(a[k][i])) {
                max = a[k][i]
                maxk = i
            }
        }
        if (maxk != k) {
            a.swapColumns(k, maxk)
            xIndexes.swap(k, maxk)
        }
        for (j in k until x.size) {
            a[k][j] /= max
        }
        b[k] /= max
        for (i in k + 1 until x.size) {
            for (j in k + 1 until x.size) {
                a[i][j] -= a[i][k] * a[k][j]
            }
            b[i] -= a[i][k] * b[k]
            a[i][k] = 0.0
        }
    }

    for (i in x.size - 1 downTo 0) {
        index = xIndexes[i].toInt()
        x[index] = b[i]
        for (j in i + 1 until x.size) {
            x[index] -= a[i][j] * x[xIndexes[j].toInt()]
        }
    }

    return x
}

fun interpolateLagrange(nodes: Vector, values: Vector, x: Double): Double {
    var result = 0.0
    var tmp: Double
    for (i in 0 until nodes.size) {
        tmp = values[i]
        for (j in 0 until nodes.size) {
            tmp *= if (i != j) (x - nodes[j]) / (nodes[i] - nodes[j]) else 1.0
        }
        result += tmp
    }
    return result
}

fun calcIntegralLeftRectangles(index: Int, nodes: Vector, solution: Vector, step: Double): Double {
    return step * (0..(nodes.size - 2)).sumByDouble { calcKernel(nodes[index], nodes[it]) * solution[it] }
}

fun calcIntegralMiddleRectangles(index: Int, nodes: Vector, solution: Vector, step: Double): Double {
    return step * (0..(nodes.size - 2)).sumByDouble { calcKernel(nodes[index], nodes[it] + step / 2) * interpolateLagrange(nodes, solution, nodes[it] + step / 2) }
}

fun calcIntegralMiddleRectanglesVolterra(index: Int, nodes: Vector, solution: Vector, step: Double) : Double {
    var solutionMiddle = Vector(solution.size)
    for (i in 0 until solution.size) {
        solutionMiddle[i] = interpolateLagrange(nodes, solution, nodes[i] + step / 2)
    }
    return calcIntegralMiddleRectangles(index, nodes.getPartVector(0, index), solutionMiddle.getPartVector(0, index), step)
}

fun calcMechanicalQuadratureFredholm(nodes: Vector, lambda: Double, intervalBottom: Double, intervalUpper: Double): Vector {
    val funcVect = Vector(nodes.size)
    val mtr = Matrix(nodes.size, nodes.size)
    val coefs = Vector(nodes.size)

    for (i in 0 until coefs.size - 1) {
        coefs[i] = (intervalUpper - intervalBottom) / nodes.size
    }
    for (i in 0 until funcVect.size) {
        funcVect[i] = calcFunction(nodes[i])
    }
    for (i in 0 until mtr.lineNum) {
        for (j in 0 until mtr.columnNum) {
            mtr[i][j] = if (i == j) 1 - lambda * coefs[j] * calcKernel(nodes[i], nodes[j]) else
                -lambda * coefs[j] * calcKernel(nodes[i], nodes[j])
        }
    }

    return gauss(mtr, funcVect)
}

fun calcMechanicalQuadratureVolterra(nodes: Vector, lambda: Double, intervalBottom: Double, intervalUpper: Double): Vector {
    val result = Vector(nodes.size)
    val step = (intervalUpper - intervalBottom) / (nodes.size - 1)
    var sum: Double

    for (i in 0 until nodes.size) {
        sum = 0.0
        for (j in 0 until i) {
            sum += calcKernel(nodes[i], nodes[j]) * step * result[j]
        }
        result[i] = (calcFunction(nodes[i]) + lambda * sum) / (1 - lambda * step * calcKernel(nodes[i], nodes[i]))
    }
    return result
}

fun calcSequentialApproximateFredholm(nodes: Vector, lambda: Double, intervalBottom: Double, intervalUpper: Double): Vector {
    var solutionCur: Vector
    val solutionNext = Vector(nodes.size)
    val step = (intervalUpper - intervalBottom) / (nodes.size - 1)
    val accuracy = Math.pow(step, 2.0)

    for (i in 0 until solutionNext.size) {
        solutionNext[i] = calcFunction(nodes[i])
    }
    do {
        solutionCur = Vector(solutionNext)
        for (i in 0 until solutionNext.size) {
            solutionNext[i] = lambda * calcIntegralLeftRectangles(i, nodes, solutionCur, step) + calcFunction(nodes[i])
        }
    } while ((solutionNext - solutionCur).norm() >= accuracy)

    return solutionNext
}

fun calcSequentialApproximateVolterra(nodes: Vector, lambda: Double, intervalBottom: Double, intervalUpper: Double): Vector {
    var solutionCur: Vector
    val solutionNext = Vector(nodes.size)
    val step = (intervalUpper - intervalBottom) / (nodes.size - 1)
    val accuracy = Math.pow(step, 2.0)

    for (i in 0 until solutionNext.size) {
        solutionNext[i] = calcFunction(nodes[i])
    }
    do {
        solutionCur = Vector(solutionNext)
        for (i in 0 until solutionNext.size) {
            solutionNext[i] = lambda * calcIntegralMiddleRectanglesVolterra(i, nodes, solutionCur, step) + calcFunction(nodes[i])
        }
    } while ((solutionNext - solutionCur).norm() >= accuracy)

    return solutionNext
}

fun main(args: Array<String>) {
    val nodeNum = 10
    val lambda = 0.5
    val intervalBottom = 0.0
    val intervalUpper = Math.PI / 2
    val step = (intervalUpper - intervalBottom) / nodeNum
    val nodes = Vector(nodeNum + 1)

    for (i in 0..nodeNum) {
        nodes[i] = intervalBottom + i * step
    }

    println("Метод механических квадратур, интегральное уравнение Фредгольма:")
    println("Узлы:")
    nodes.print()
    println("Значения искомой функции в узлax:")
    calcMechanicalQuadratureFredholm(nodes, lambda, intervalBottom, intervalUpper).print()
    println()

    println("Метод механических квадратур, интегральное уравнение Вольтерра:")
    println("Узлы:")
    nodes.print()
    println("Значения искомой функции в узлax:")
    calcMechanicalQuadratureVolterra(nodes, lambda, intervalBottom, intervalUpper).print()
    println()

    println("Метод последовательных приближений, интегральное уравнение Фредгольма:")
    println("Узлы:")
    nodes.print()
    println("Значения искомой функции в узлax:")
    calcSequentialApproximateFredholm(nodes, lambda, intervalBottom, intervalUpper).print()
    println()

    println("Метод последовательных приближений, интегральное уравнение Вольтерра:")
    println("Узлы:")
    nodes.print()
    println("Значения искомой функции в узлax:")
    calcSequentialApproximateVolterra(nodes, lambda, intervalBottom, intervalUpper).print()
    println()
}