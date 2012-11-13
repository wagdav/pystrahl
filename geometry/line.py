"""
Bresenham's line drawing algorithm.  The implementation is based on the sample
given at [1].

[1] http://rosettacode.org/wiki/Bitmap/Bresenham's_line_algorithm

>>> bitmap = Bitmap(17,17)
>>> for points in ((1,8,8,16),(8,16,16,8),(16,8,8,1),(8,1,1,8)):
...     bitmap.line(*points)
>>> bitmap.chardisplay()
+-----------------+
|        @        |
|       @ @       |
|      @   @      |
|     @     @     |
|    @       @    |
|    @        @   |
|   @          @  |
|  @            @ |
| @              @|
|  @            @ |
|   @          @  |
|    @       @@   |
|     @     @     |
|      @   @      |
|       @ @       |
|        @        |
|                 |
+-----------------+
"""

def line(x0, y0, x1, y1):
    line_points = []
    dx = abs(x1 - x0)
    dy = abs(y1 - y0)
    x, y = x0, y0
    sx = -1 if x0 > x1 else 1
    sy = -1 if y0 > y1 else 1
    if dx > dy:
        err = dx / 2.0
        while x != x1:
            line_points.append((x,y))
            err -= dy
            if err < 0:
                y += sy
                err += dx
            x += sx
    else:
        err = dy / 2.0
        while y != y1:
            line_points.append((x,y))
            err -= dx
            if err < 0:
                x += sx
                err += dy
            y += sy
    line_points.append((x,y))
    return line_points


class Bitmap():
    def __init__(self, width=40, height=40):
        assert width > 0 and height > 0
        self.width = width
        self.height = height
        self.map = set()

    def chardisplay(self):
        aa = [[ (w,h) for w in xrange(self.width)] for h in xrange(self.height)]
        txt = [''.join('@' if (bit in self.map) else ' ' for bit in row) for row in aa]
        # Boxing
        txt = ['|'+row+'|' for row in txt]
        txt.insert(0, '+' + '-' * self.width + '+')
        txt.append('+' + '-' * self.width + '+')
        print('\n'.join(reversed(txt)))

    def set(self, x, y):
        self.map.add((x, y))

    def get(self, x, y):
        return (x, y) in self.map

    def line(self, x0, y0, x1, y1):
        self.map.update(line(x0, y0, x1, y1))


if __name__ == '__main__':
    import doctest
    doctest.testmod()
