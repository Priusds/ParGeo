import math


def normalize(point):
    norm = math.sqrt(sum([point[i] ** 2 for i in range(3)]))
    return tuple([point[i] / norm for i in range(3)])


def locate(point, normalVector, supportVector):
    """
    True if in else False
    """
    s = sum([(point[i] - supportVector[i]) * normalVector[i] for i in range(3)])
    return s <= 0


def merge_loops(loop0, loop1, h0, h1, geo):
    # loop = (id , [lines])

    # assert len(intersectionPoints) == 2:
    # loop0_clipped = [] #[line for line in loop0[1] if line in geo["lines"]]
    loop = []
    loop1_added = False
    start_loop1 = None
    for i, line in enumerate(loop0[1]):
        if line in geo["lines"]:
            loop.append(line)

        elif not loop1_added:
            start_loop1 = len(loop)
            if len(loop) > 0:
                start0, end0 = geo["lines"][abs(loop[-1])]
                last_vidx = start0 if loop[-1] < 0 else end0
                if last_vidx in geo["lines"][abs(loop1[1][0])]:
                    loop += loop1[1]
                else:
                    loop += [-l for l in reversed(loop1[1])]
                loop1_added = True
    end_loop1 = start_loop1 + len(loop1)
    return (
        loop,
        loop[start_loop1 : end_loop1 + 2],
        loop[end_loop1 + 2 :] + loop[: start_loop1 - 1],
    )
