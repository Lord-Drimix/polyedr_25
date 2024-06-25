from math import pi
from functools import reduce
from operator import add
from common.r3 import R3
import math


class Segment:
    """ Одномерный отрезок """
    # Параметры конструктора: начало и конец отрезка (числа)

    def __init__(self, beg, fin):
        self.beg, self.fin = beg, fin

    # Отрезок вырожден?
    def is_degenerate(self):
        return self.beg >= self.fin

    # Пересечение с отрезком
    def intersect(self, other):
        if other.beg > self.beg:
            self.beg = other.beg
        if other.fin < self.fin:
            self.fin = other.fin
        return self

    # Разность отрезков
    # Разность двух отрезков всегда является списком из двух отрезков!
    def subtraction(self, other):
        return [Segment(
            self.beg, self.fin if self.fin < other.beg else other.beg),
            Segment(self.beg if self.beg > other.fin else other.fin, self.fin)]


class Edge:
    """ Ребро полиэдра """
    # Начало и конец стандартного одномерного отрезка
    SBEG, SFIN = 0.0, 1.0

    # Параметры конструктора: начало и конец ребра (точки в R3)
    def __init__(self, beg, fin):
        self.beg, self.fin = beg, fin
        # Список «просветов»
        self.gaps = [Segment(Edge.SBEG, Edge.SFIN)]
        self.shaded = False

    # Учёт тени от одной грани
    def shadow(self, facet):
        # «Вертикальная» грань не затеняет ничего
        if facet.is_vertical():
            return
        # Нахождение одномерной тени на ребре
        shade = Segment(Edge.SBEG, Edge.SFIN)
        for u, v in zip(facet.vertexes, facet.v_normals()):
            shade.intersect(self.intersect_edge_with_normal(u, v))#находим пересечение с вертикальными полупространствами(во всем форе,но shade intersect непосредственно)
            #будут пересекаться тени граней
            if shade.is_degenerate():
                return

        shade.intersect(
            self.intersect_edge_with_normal(
                facet.vertexes[0], facet.h_normal()))#для горизонтальных граней
        if shade.is_degenerate():
            return
        # Преобразование списка «просветов», если тень невырождена
        gaps = [s.subtraction(shade) for s in self.gaps] #вычитание теней из просвета
        self.gaps = [
            s for s in reduce(add, gaps, []) if not s.is_degenerate()]# двумерный в одномерный и убираем вырожденные тени

    # Преобразование одномерных координат в трёхмерные
    def r3(self, t):
        return self.beg * (Edge.SFIN - t) + self.fin * t

    # Пересечение ребра с полупространством, задаваемым точкой (a)
    # на плоскости и вектором внешней нормали (n) к ней
    def intersect_edge_with_normal(self, a, n): #считаем скалярное произведение  Находит пересечение вектора с (в какой точке плоскость делит данный отрезок)
        f0, f1 = n.dot(self.beg - a), n.dot(self.fin - a) # f1 и f0 числа, знак которых показывает с какой стороны от плоскости лежит ребро и анализируем знаки
        if f0 >= 0.0 and f1 >= 0.0: # делаем ребро вырожденнным, то мы меняем ребро начало и конец. Значит все ребро затемнено
            return Segment(Edge.SFIN, Edge.SBEG)
        if f0 < 0.0 and f1 < 0.0:
            return Segment(Edge.SBEG, Edge.SFIN)
        x = - f0 / (f1 - f0) #точка,  которой плоскость пересекает ребро
        return Segment(Edge.SBEG, x) if f0 < 0.0 else Segment(x, Edge.SFIN)


class Facet:
    """ Грань полиэдра """
    # Параметры конструктора: список вершин

    def __init__(self, vertexes, edges=[]):
        self.vertexes = vertexes
        self.edges = edges

    # «Вертикальна» ли грань?
    def is_vertical(self):
        return self.h_normal().dot(Polyedr.V) == 0.0

    # Нормаль к «горизонтальному» полупространству
    def h_normal(self):
        n = (
            self.vertexes[1] - self.vertexes[0]).cross(
            self.vertexes[2] - self.vertexes[0])
        return n * (-1.0) if n.dot(Polyedr.V) < 0.0 else n

    # Нормали к «вертикальным» полупространствам, причём k-я из них
    # является нормалью к грани, которая содержит ребро, соединяющее
    # вершины с индексами k-1 и k
    def v_normals(self):
        return [self._vert(x) for x in range(len(self.vertexes))]

    # Вспомогательный метод
    def _vert(self, k):
        n = (self.vertexes[k] - self.vertexes[k - 1]).cross(Polyedr.V)
        return n * \
            (-1.0) if n.dot(self.vertexes[k - 1] - self.center()) < 0.0 else n
    #непосредственно считает нормали к вертикальным

    # Центр грани
    def center(self):
        return sum(self.vertexes, R3(0.0, 0.0, 0.0)) * \
            (1.0 / len(self.vertexes))


class Polyedr:
    """ Полиэдр """
    # вектор проектирования
    V = R3(0.0, 0.0, 1.0)

    # Параметры конструктора: файл, задающий полиэдр
    def __init__(self, file):

        # списки вершин, рёбер и граней полиэдра
        self.vertexes, self.edges, self.facets = [], [], []
        self.perimeter = 0

        # список строк файла
        with open(file) as f:
            for i, line in enumerate(f):
                if i == 0:
                    # обрабатываем первую строку; buf - вспомогательный массив
                    buf = line.split()
                    # коэффициент гомотетии
                    self.c = float(buf.pop(0))
                    # углы Эйлера, определяющие вращение
                    alpha, beta, gamma = (float(x) * pi / 180.0 for x in buf)
                elif i == 1:
                    # во второй строке число вершин, граней и рёбер полиэдра
                    nv, nf, ne = (int(x) for x in line.split())
                elif i < nv + 2:
                    # задание всех вершин полиэдра
                    x, y, z = (float(x) for x in line.split())
                    self.vertexes.append(R3(x, y, z).rz(
                        alpha).ry(beta).rz(gamma) * self.c)
                else:
                    # вспомогательный массив
                    buf = line.split()
                    # количество вершин очередной грани
                    size = int(buf.pop(0))
                    # массив вершин этой грани
                    vertexes = list(self.vertexes[int(n) - 1] for n in buf)
                    edges = []
                    # задание рёбер грани
                    for n in range(size):
                        edge = Edge(vertexes[n - 1], vertexes[n])
                        edges.append(edge)
                        self.edges.append(edge)
                    # задание самой грани
                    self.facets.append(Facet(vertexes, edges))

    def shad(self): #сичтаем тени на каждом ребре первые два фора
        for e in self.edges:
            for f in self.facets:
                e.shadow(f)

            sigma = 0 #и ниже наъходим полностью затемненные ребра, на них либо нет просветов либо Сигма
            for x in e.gaps: #по всем просвтам, фактически считаю длину каждого просвета
                sigma += x.fin - x.beg
            if len(e.gaps) == 0 or sigma < 1e-4:
                e.shaded = True

        for k in self.facets: #фор по всем граням
            shaded_facet = True #предполагаем, что гшрань полностью затемнена
            for e in k.edges: #фором проверяем, что все ее ребра затемнены.
                shaded_facet = shaded_facet and e.shaded #тру если грань полностью затемнена,
            if shaded_facet:
                if math.sqrt(k.center().x ** 2 + k.center().y ** 2
                             + k.center().z ** 2) < 2 * self.c:
                    for e in k.edges: #добавляем длину проекции
                        self.perimeter += math.sqrt((e.beg.x - e.fin.x) ** 2 +
                                                    (e.beg.y - e.fin.y) ** 2)
        self.perimeter = self.perimeter / self.c
        return self.perimeter

    # Метод изображения полиэдра
    def draw(self, tk):  # pragma: no cover
        tk.clean()
        for e in self.edges:
            for f in self.facets:
                e.shadow(f)
            for s in e.gaps:
                tk.draw_line(e.r3(s.beg), e.r3(s.fin))

            sigma = 0
            for x in e.gaps:
                sigma += x.fin - x.beg
            if len(e.gaps) == 0 or sigma < 1e-4:
                e.shaded = True

        for k in self.facets:
            shaded_facet = True
            for e in k.edges:
                shaded_facet = shaded_facet and e.shaded
            if shaded_facet:
                if math.sqrt(k.center().x**2 + k.center().y**2 +
                             k.center().z**2) < 2 * self.c:
                    for e in k.edges:
                        self.perimeter += math.sqrt((e.beg.x - e.fin.x) ** 2 +
                                                    (e.beg.y - e.fin.y) ** 2)
        self.perimeter = self.perimeter / self.c
