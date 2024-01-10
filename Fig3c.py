import cv2
import numpy as np
import os

DIR = "Data/Fig3c/"
Files = os.listdir(DIR)
Files = sorted([item for item in Files if "Mock" not in item and item != '.DS_Store'])

output_folder = "Output_slim"
if not os.path.exists(output_folder):
    os.mkdir(output_folder)

for file in Files:
    image_read = cv2.imread(os.path.join(DIR, file))

    image_read_calc = np.copy(image_read)
    image_read_calc = cv2.GaussianBlur(image_read_calc, (9, 9), cv2.BORDER_DEFAULT)
    image_read_calc[image_read_calc[:, :, 0] / 1.1 < np.max(image_read_calc[:, :, 1:3], axis=-1)] = 0
    image_read_calc[image_read_calc[:, :, 0] < 130] = 0

    imgray = cv2.cvtColor(image_read_calc, cv2.COLOR_BGR2GRAY)
    block_size = int(round(len(imgray) / 3))
    if block_size % 2 == 0:
        block_size += 1
    block_size = int(block_size)
    imgray = cv2.adaptiveThreshold(imgray, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY, block_size, -2)

    contours, hierarchy = cv2.findContours(imgray, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
    contours = sorted(contours, key=cv2.contourArea, reverse=True)
    contours = [contour for contour in contours if cv2.contourArea(contour) > 100]

    m0 = cv2.moments(contours[0])
    cx0 = int(m0["m10"] / m0["m00"])
    cy0 = int(m0["m01"] / m0["m00"])
    coord0 = np.array([cx0,cy0])
    combine = [0]
    for index, contour in enumerate(contours):
        m = cv2.moments(contour)
        cx = int(m["m10"] / m["m00"])
        cy = int(m["m01"] / m["m00"])
        coord = np.array([cx, cy])
        if np.linalg.norm(coord0 - coord) < 500:
            combine.append(index)

    contours = [contour for index, contour in enumerate(contours) if index in combine]
    contours = np.array([np.vstack(contours)])

    rect = cv2.minAreaRect(contours[0])
    box = cv2.boxPoints(rect)
    box = np.int0(box)
    y = [item for item in box[:, 1] if item > 0]
    x = box[:, 0]

    image_read = image_read[np.min(y):np.max(y), np.min(x):np.max(x)]

    image_read_calc = np.copy(image_read)
    image_read_calc[image_read_calc[:, :, 0] - 10 < np.max(image_read_calc[:, :, 1:3], axis=-1)] = 0
    image_read_calc[np.sum(image_read_calc[:, :, 0:3], axis=-1) > 450] = 0

    area2 = cv2.countNonZero(cv2.cvtColor(image_read_calc, cv2.COLOR_BGR2GRAY))
    cv2.imwrite(f"{output_folder}/{file}.png", image_read_calc)

    print(f"{file},{area2}")