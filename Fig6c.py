import numpy as np
import math
import os
import cv2
import scipy
from scipy import ndimage, signal
import pandas as pd
from operator import itemgetter

dir = "Data/Fig6"
files = os.listdir(dir)
files = sorted([item for item in files if item.endswith(".tif")])

result_path = "Fig6_Output"
if not os.path.exists(result_path):
    os.makedirs(result_path)

def remove_neighbours(distance, maximum_coords_y, maximum_coords_x, image_dist):
    index = range(len(maximum_coords_x))
    values = image_dist[(maximum_coords_y, maximum_coords_x)]

    maximum_coords_y_rounded = np.rint((np.rint(np.array(maximum_coords_y) / distance) * distance)).astype(np.uint)
    maximum_coords_x_rounded = np.rint((np.rint(np.array(maximum_coords_x) / distance) * distance)).astype(np.uint)
    maximum_coords_x_rounded[maximum_coords_x_rounded >= len(image_dist[0])] = len(image_dist[0]) - 1
    maximum_coords_y_rounded[maximum_coords_y_rounded >= len(image_dist)] = len(image_dist) - 1
    rounded_list = sorted(zip(maximum_coords_y_rounded, maximum_coords_x_rounded, values, index), key=itemgetter(2),
                          reverse=True)

    seen = set()
    seen2 = set()
    [seen.add(x[:2]) or seen2.add(x[3]) for x in rounded_list if x[:2] not in seen]
    maximum_coords_x = [item for index, item in enumerate(maximum_coords_x) if index in seen2]
    maximum_coords_y = [item for index, item in enumerate(maximum_coords_y) if index in seen2]
    return maximum_coords_x, maximum_coords_y

for file_name in files:
    file_name_output = os.path.join(result_path,file_name)
    file_name = os.path.join(dir,file_name)

    im = cv2.imread(file_name, 0)
    im = cv2.GaussianBlur(im, (3, 3), 0)
    number = 2
    im[im < number] = 0
    im[im >= number] = 255

    x_data = np.max(im, axis=0)
    y_data = np.max(im, axis=1)
    x_data_min = np.argmax(x_data)
    x_data_max = len(x_data) - np.argmax(x_data[::-1])
    y_data_min = np.argmax(y_data)
    y_data_max = len(y_data) - np.argmax(y_data[::-1])

    thickness = 30
    side = False
    top = False

    if x_data_min < thickness:
        side = "Left"
        im = im[:, x_data_min:]
    elif x_data_max > len(im[0]) - thickness:
        side = "Right"
        im = im[:, :x_data_max]
    if y_data_min < thickness:
        top = "top"
        im = im[y_data_min:, :]
    elif y_data_min > len(im) - thickness:
        top = "bottom"
        im = im[:y_data_min, :]

    print(file_name,side,top)

    if side and top:
        #Draw boundries x-axis corners
        if len(np.where(im[0] == 255)[0]) > 0:
            x = np.where(im[0] == 255)[0]
            if min(x) < (len(im[0]) - max(x)):
                cv2.line(im, (0, 0), (max(x), 0), 255, thickness)
            else:
                cv2.line(im, (len(im[0]), 0), (min(x), 0), 255, thickness)

        if len(np.where(im[-1] == 255)[0]) > 0:
            x = np.where(im[-1] == 255)[0]
            if min(x) < (len(im[0]) - max(x)):
                cv2.line(im, (0, len(im)), (max(x), len(im)), 255, thickness)
            else:
                cv2.line(im, (len(im[0]), len(im)), (min(x), len(im)), 255, thickness)

        #y-axis
        if len(np.where(im[:,0] == 255)[0]) > 0:
            x = np.where(im[:,0] == 255)[0]
            if min(x) < (len(im[:,0]) - max(x)):
                cv2.line(im, (0,0), (0,max(x)), 255, thickness)
            else:
                cv2.line(im, (0,len(im)), (0,min(x)), 255, thickness)

        if len(np.where(im[:,-1] == 255)[0]) > 0:
            x = np.where(im[:,-1] == 255)[0]
            if min(x) < (len(im[:,0]) - max(x)):
                cv2.line(im, (len(im[0]),0), (len(im[0]),max(x)), 255, thickness)
            else:
                cv2.line(im, (len(im[0]),len(im)), (len(im[0]),min(x)), 255, thickness)

    elif side and not top:
        #y-axis
        if len(np.where(im[:,0] == 255)[0]) > len(np.where(im[:,-1] == 255)[0]):
            x = np.where(im[:,0:5] == 255)[0]
            cv2.line(im, (0,min(x)), (0,max(x)), 255, thickness)
        else:
            x = np.where(im[:,-5:-1] == 255)[0]
            cv2.line(im, (len(im[0]),min(x)), (len(im[0]),max(x)), 255, thickness)

    im = im[y_data_min:y_data_max, x_data_min:x_data_max]

    imgray = im.copy()
    block_size = int(round(len(imgray)/3))
    if block_size % 2 == 0:
        block_size += 1
    block_size = int(block_size)
    imgray = cv2.adaptiveThreshold(imgray, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY, block_size, -2)

    contours, hierarchy = cv2.findContours(imgray, cv2.RETR_LIST, cv2.CHAIN_APPROX_SIMPLE)
    contours = sorted(contours, key=cv2.contourArea, reverse=True)
    im = np.zeros_like(im)

    cv2.drawContours(im,contours, 0, 255, -1)
    DistanceMap = scipy.ndimage.distance_transform_edt(im)

    gaussian_SharpMap = cv2.GaussianBlur(DistanceMap, (0, 0), 2.0)
    DistanceMap = cv2.addWeighted(DistanceMap, 2.0, gaussian_SharpMap, -1.0, 0)

    MaxFilter = ndimage.maximum_filter(DistanceMap, size=5, mode='constant')
    diff = MaxFilter - DistanceMap
    Mask = (diff.astype(int) == 0)
    black = np.zeros_like(im)
    black[Mask] = im[Mask]
    RAW_Root = np.where((black != 0) & (DistanceMap != 0))

    Selection = cv2.blur(black, (5,5))
    Selection[Selection > 0] = 255
    contours, = cv2.findContours(Selection, cv2.RETR_LIST, cv2.CHAIN_APPROX_SIMPLE)
    contours = sorted(contours, key=cv2.contourArea, reverse=True)#[0]

    RAW_Root_X = RAW_Root[1]
    RAW_Root_Y = RAW_Root[0]

    Selection = np.zeros_like(im)
    cv2.drawContours(Selection,contours, 0, 255, -1)
    cv2.drawContours(Selection,contours, 1, 255, -1)
    points = Selection[RAW_Root_Y,RAW_Root_X]
    RAW_Root_Y = list(RAW_Root_Y[points == 255])
    RAW_Root_X = list(RAW_Root_X[points == 255])

    #maximum_coords = [item for index, item in enumerate(maximum_coords) if cv2.pointPolygonTest(contours, (maximum_coords_y[index],maximum_coords_x[index]), True) >= 0]
    if side == "Left":
        First_point_index = np.argmax(RAW_Root_X)
    else:
        First_point_index = np.argmin(RAW_Root_X)

    First_point = RAW_Root_Y[First_point_index],RAW_Root_X[First_point_index]
    RAW_Root_X, RAW_Root_Y = remove_neighbours(5, RAW_Root_Y, RAW_Root_X, DistanceMap)

    CheckPoints = list(zip(RAW_Root_Y, RAW_Root_X))
    RAW_Root2 = [First_point]
    while len(CheckPoints) > 0:
        last = RAW_Root2[-1]
        min_dist = [math.sqrt((coordinate[1] - last[1]) ** 2 + (coordinate[0] - last[0]) ** 2) for coordinate in CheckPoints]
        new = CheckPoints[np.argmin(min_dist)]
        parts = int(round(min(min_dist)/5))
        intermediate = list(zip(np.linspace(last[0], new[0], parts + 1, dtype=np.uint64), np.linspace(last[1], new[1], parts + 1,dtype=np.uint64)))
        RAW_Root2 += intermediate
        CheckPoints.remove(new)

    seen = set()
    RAW_Root2 = [x for x in RAW_Root2 if not (tuple(x) in seen or seen.add(tuple(x)))]
    RAW_Root2 = [item for item in RAW_Root2 if DistanceMap[item[0],item[1]] > 0]

    RAW_Root2 = RAW_Root2[::-1]
    Root_Distance = [0] + [math.sqrt((int(RAW_Root2[i][1]) - int(RAW_Root2[i+1][1])) ** 2 + (int(RAW_Root2[i][0]) - int(RAW_Root2[i+1][0])) ** 2) for i in range(len(RAW_Root2)-1)]
    Root_Distance = np.cumsum(Root_Distance)

    x = [item[1] for item in RAW_Root2]
    y = [item[0] for item in RAW_Root2]

    Width = [DistanceMap[item[0],item[1]] for item in RAW_Root2]
    dict = {"distance": Root_Distance, "length": Root_Distance, "Width": Width, "x": x, "y": y}
    df = pd.DataFrame.from_dict(dict)
    df.to_csv(f"{file_name_output}.csv")

    DistanceMap = DistanceMap*255/np.max(DistanceMap)
    DistanceMap_RGB = np.stack((DistanceMap,) * 3, axis=-1)
    colors = np.linspace((255, 0, 0),(0, 0, 255), len(RAW_Root2))

    for index, item in enumerate(RAW_Root2):
        cv2.circle(DistanceMap_RGB, (int(item[1]), int(item[0])), 3, color=colors[index], thickness=-1,)

    cv2.imwrite(f"{file_name_output}_DistanceMap.png", DistanceMap_RGB)