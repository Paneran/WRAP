# pyserial validation
import serial
import csv

ser = serial.Serial('COM4', 115200, timeout=30, parity=serial.PARITY_EVEN, rtscts=1)
s = ser.read(10000)       # read up to one hundred bytes
                          # or as much is in the buffer

f = open('./data/stm_costasamples.csv', 'w')

# create the csv writer
writer = csv.writer(f)

# write a row to the csv file
writer.writerow(s)

# close the file
f.close()