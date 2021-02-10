import tkinter as tk
from tkinter import *
from tkinter import filedialog
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.io import fits
from tkinter import simpledialog
from tkinter import messagebox as msb
import numpy as np
import astropy.units as u
from astropy.wcs import WCS
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clipped_stats
from photutils import DAOStarFinder
from astroquery.vizier import Vizier
from astropy.nddata import Cutout2D
import math
import time
from matplotlib.colors import LogNorm
from matplotlib.widgets import Slider, Button
import csv
root= tk.Tk()
root.title('FITS')
root.geometry("800x800+50+150")
img=PhotoImage(file='instruction.png')
pixels_over_limit = 0
pixels_over_limit_in_reference_stars = 0
FITS_counter = 0
stars_coordinates = []
main_star_brightness = []
main_background_brightness = []
pomiar_gwiazdy_znormalizowany = []
main_star_brigtness_1_debug = []
reference_ring_values = []
first_fits_references_stars = []
reference_stars_limit = 0
srodek_x = .0
srodek_y = .0
button = 0
first_ref = 0

csv_filenames = []
csv_main_values = []
csv_background_values = []
csv_main_brightness_values = []
csv_ref_coords = []
csv_ring_values = []
csv_ref_ring_values = []

def set_srodek_y(var):
    global srodek_y
    srodek_y = int(var)

def set_srodek_x(var):
    global srodek_x
    srodek_x = int(var)

def get_cursor_position(event):
    cursor_x,cursor_y,button_x = event.xdata, event.ydata, event.button
    if(button_x == 3):
        set_srodek_x(cursor_x)
        set_srodek_y(cursor_y)

def open_new_window():
    R=0
    newWindow = Toplevel(root)
    newWindow.title("Instruction")
    newWindow.geometry("799x695+800+200")
    Label(newWindow,image=img).pack()
    r = int(simpledialog.askstring(title="r", prompt="Podaj r:"))
    while R<=(2*r)+5:
        R = int(simpledialog.askstring(title="R", prompt="Podaj R:"))
    if r<=R: newWindow.destroy()
    return r,R

def burned_out_pixel():
    print("Znaleziono przepalony piksel, ustalenie tranzytu jest niemożliwe")

def ring_estimation(file,r,R,srodek_x,srodek_y,main_counter):
    file = file
    hdu = fits.open(file)
    scidata = hdu[0].data
    scidata_height = scidata.shape[0]
    scidata_width = scidata.shape[1]
    hud = fits.open(file)[0]
    im_data = hud.data
    scidata = im_data
    limit = 65000 #limit przepalenia
    matryca = []
    srodek = []
    drawing=np.empty((scidata_height,scidata_width))
    drawing2=np.empty((scidata_height,scidata_width))
    size = R
    for i in range(0,size+1):
        for k in range(0,size+1):
            x = srodek_x + i
            y = srodek_y + k
            if (math.sqrt(i*i+k*k))>2*r and ((math.sqrt(i*i+k*k)))<R:
                    drawing[x][y] = scidata[x,y]
                    matryca.append(scidata[x,y])
        for l in range(0,size+1):
            x = srodek_x + i
            y = srodek_y - l
            if (math.sqrt(i*i+l*l))>2*r and (math.sqrt(i*i+l*l))<R:
                    matryca.append(im_data[x,y])
                    drawing[x][y] = scidata[x,y]
        for k in range(0,size+1):
            x = srodek_x - i
            y = srodek_y + k
            if (math.sqrt(i*i+k*k))>2*r and (math.sqrt(i*i+k*k))<R:
                    matryca.append(im_data[x,y])
                    drawing[x][y] = scidata[x,y]
        for l in range(0,size+1):
            x = srodek_x - i
            y = srodek_y - l
            if (math.sqrt(i*i+l*l))>2*r and (math.sqrt(i*i+l*l))<R:
                    matryca.append(im_data[x,y])
                    drawing[x][y] = scidata[x,y]
    burned_pixel_status = False
    drawing2=np.empty((scidata_height,scidata_width))
    for i in range(0,r+1):
        for k in range(0,r+1):
            x = srodek_x + i
            y = srodek_y + k
            if (math.sqrt(i*i+k*k))<r:
                if(scidata[x,y]<limit):
                    srodek.append(im_data[x,y])
                    drawing2[x][y] = scidata[x,y]
                else:
                    burned_out_pixel()
                    burned_pixel_status = True
                    break
        if(burned_pixel_status==True): break
        for l in range(0,r+1):
            x = srodek_x + i
            y = srodek_y - l
            if (math.sqrt(i*i+l*l))<r:
                if(scidata[x,y]<limit):
                    srodek.append(im_data[x,y])
                    drawing2[x][y] = scidata[x,y]
                else:
                    burned_out_pixel()
                    burned_pixel_status = True
                    break
        if(burned_pixel_status==True): break
        for k in range(0,r+1):
            x = srodek_x - i
            y = srodek_y + k
            if (math.sqrt(i*i+k*k))<r:
                if(scidata[x,y]<limit):
                    srodek.append(im_data[x,y])
                    drawing2[x][y] = scidata[x,y]
                else:
                    burned_out_pixel()
                    burned_pixel_status = True
                    break
        if(burned_pixel_status==True): break
        for l in range(0,r+1):
            x = srodek_x - i
            y = srodek_y - l
            if (math.sqrt(i*i+l*l))<r:
                if(scidata[x,y]<limit):
                    srodek.append(im_data[x,y])
                    drawing2[x][y] = scidata[x,y]
                else:
                    burned_out_pixel()
                    burned_pixel_status = True
                    break
        if(burned_pixel_status==True): break
    if(burned_pixel_status == False):
        #Tymczasowe
        gwiazda = sum(srodek)
        tlo = len(srodek) * np.mean(matryca)
        wartosc_gwiazdy = gwiazda - tlo
        #testuje to
        if main_counter == 1:pass
            # plt.imshow(drawing)
            # plt.show()
            # plt.imshow(drawing2)
            # plt.show()
        hdu.close()
        return (wartosc_gwiazdy/100,tlo)

def reference_ring(file,coord_x,coord_y,tlo,main_counter):
    hdu = fits.open(file)
    scidata = hdu[0].data
    scidata_height = scidata.shape[0]
    scidata_width = scidata.shape[1]
    im_data = hdu[0].data
    r = 15 #stały promień dla gwiazd odniesienia
    limit = 65000 #limit przepalenia
    srodek = []
    drawing2=np.empty((scidata_height,scidata_width))
    burned_pixel_status = False
    srodek_x = coord_y
    srodek_y = coord_x
    for i in range(0,r+1):
        for k in range(0,r+1):
            x = srodek_x + i
            y = srodek_y + k
            if (math.sqrt(i*i+k*k))<r:
                if(x<scidata_width and y<scidata_height):
                    if(scidata[x,y]<limit):
                        srodek.append(im_data[x,y])
                        drawing2[x][y] = scidata[x,y]
                    else:
                        burned_out_pixel()
                        burned_pixel_status = True
                        break
        if(burned_pixel_status==True): break
        for l in range(0,r+1):
            x = srodek_x + i
            y = srodek_y - l
            if (math.sqrt(i*i+l*l))<r:
                if(x<scidata_width and y<scidata_height):
                    if(scidata[x,y]<limit):
                        srodek.append(im_data[x,y])
                        drawing2[x][y] = scidata[x,y]
                    else:
                        burned_out_pixel()
                        burned_pixel_status = True
                        break
        if(burned_pixel_status==True): break
        for k in range(0,r+1):
            x = srodek_x - i
            y = srodek_y + k
            if (math.sqrt(i*i+k*k))<r:
                if(x<scidata_width and y<scidata_height):
                    if(scidata[x,y]<limit):
                        srodek.append(im_data[x,y])
                        drawing2[x][y] = scidata[x,y]
                    else:
                        burned_out_pixel()
                        burned_pixel_status = True
                        break
        if(burned_pixel_status==True): break
        for l in range(0,r+1):
            x = srodek_x - i
            y = srodek_y - l
            if (math.sqrt(i*i+l*l))<r:
                if(x<scidata_width and y<scidata_height):
                    if(scidata[x,y]<limit):
                        srodek.append(im_data[x,y])
                        drawing2[x][y] = scidata[x,y]
                    else:
                        burned_out_pixel()
                        burned_pixel_status = True
                        break
        if(burned_pixel_status==True):
            break
            msb.showinfo("Error","Znaleziono przepalone piksele, określenie tranzytu jest niemożliwe")
    if(burned_pixel_status == False):
        gwiazda = sum(srodek)
        wartosc_gwiazdy = gwiazda-tlo
        hdu.close()
        drawing2
        return (gwiazda)

def reference_star_finder (file,x,y):
    pixel_stars = []
    master_fits = fits.open(file)[0]
    wcs = WCS(master_fits.header)
    data = master_fits.data
    position = (x, y)
    srodek_x = x
    srodek_y = y
    size_x = (data.shape[1]-srodek_x)
    size_y = (data.shape[0]-srodek_y)
    size = (size_x, size_y)

    cutout = Cutout2D(master_fits.data, position=position, size=size, wcs=wcs)
    temp_fits = cutout.data
    mean, median, std = sigma_clipped_stats(temp_fits)
    star_finder = DAOStarFinder(threshold=8*std, fwhm=3)
    pixel_stars = star_finder(temp_fits - median)
    print("ilość gwiazd odniesienia przed poprawieniem: ",len(pixel_stars))
    # if(data.shape[1]-srodek_x>700 and data.shape[0]-srodek_y>700):
    #     size = (700,700)
    #     cutout2 = Cutout2D(master_fits.data, position=position, size=size, wcs=wcs)
    #     temp_fits2 = cutout2.data
    #     mean, median, std = sigma_clipped_stats(temp_fits2)
    #     star_finder2 = DAOStarFinder(threshold=8*std, fwhm=4)
    #     pixel_stars2 = star_finder2(data - median)
    #     print("ilość gwiazd odniesienia po poprawieniu: ",len(pixel_stars2))
    #     return(pixel_stars2)
    # else: return(pixel_stars)
    return(pixel_stars)

def open_fits():
    fileName = filedialog.askopenfilename(title = "Select file",
    filetypes = (("FIT files","*.fit"),("FITS files","*.fits")))
    if fileName:
        file = fits.open(fileName)[0]
        img_data = file.data
        fig = plt.figure()
        plt.imshow(img_data,cmap='gray')
        plt.colorbar()
        fig.suptitle(fileName)
        canvas = FigureCanvasTkAgg(fig, master=root)
        toolbar = NavigationToolbar2Tk(canvas, root, pack_toolbar=False)
        toolbar.pack(side=tk.BOTTOM, fill=tk.X)
        toolbar.update()
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

def calibration(file,main_counter):
    opened_file = fits.open(file)
    msb.showinfo("Informacja","Teraz podaj ścieżkę dostępu do klatki ciemnej")
    dark = filedialog.askopenfilename(title = "Select file",
    filetypes = (("FITS files","*.fits"),("FIT files","*.fit")))
    msb.showinfo("Informacja","Teraz podaj ścieżkę dostępu do klatki płaskiej")
    flat = filedialog.askopenfilename(title = "Select file",
    filetypes = (("FITS files","*.fits"),("FIT files","*.fit")))
    wcs = WCS(opened_file[0].header)
    dark_temp = fits.open(dark)
    flat_temp = fits.open(flat)
    data = opened_file[0].data
    dark_data = dark_temp[0].data
    flat_data = flat_temp[0].data
    print("data: ",data)
    print("dark data:",dark_data)
    print("flat data:",flat_data)
    sub_FITS = np.subtract(data, dark_data)
    calib_FITS = np.divide(sub_FITS, flat_data)
    CalibFITS = fits.PrimaryHDU(calib_FITS)
    calibrated_filename = '../data/calibratedFITS/CalibFITS_'+str(main_counter)+'.fit'
    hdr = opened_file[0].header
    hdu = CalibFITS.data
    hdr.update(wcs.to_header())
    fits.writeto(calibrated_filename,hdu,hdr, overwrite=True)
    if main_counter == 0:
        file_show = fits.open(file)[0]
        img_data = file_show.data #to jest raw
        c_min = 2000
        c_max = 4500
        axcolor = 'lightgoldenrodyellow'
        after_calib_data = hdu
        plot, ax = plt.subplots(1, 2, figsize = (15, 7))
        plot.canvas.set_window_title('Kalibracja')
        img = ax[0].imshow(img_data, cmap='gray')
        ax[0].set_title("BEFORE")
        divider = make_axes_locatable(ax[0])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plot.colorbar(img, cax = cax)
        img2 = ax[1].imshow(after_calib_data, cmap='gray')
        divider2 = make_axes_locatable(ax[1])
        cax2 = divider2.append_axes("right", size="5%", pad=0.05)
        ax[1].set_title("AFTER")
        plot.colorbar(img2, cax=cax2)
        ax_cmin = plt.axes([0.05, 0.05, 0.7, 0.03]) #left,bottom,width,height
        ax_cmax  = plt.axes([0.05, 0.02, 0.7, 0.03])
        s_cmin = Slider(ax_cmin, 'white', 0, 65000, valinit=c_min)
        s_cmax = Slider(ax_cmax, 'black', 0, 65000, valinit=c_max)
        def update(val, s=None):
            _cmin = s_cmin.val
            _cmax = s_cmax.val
            img.set_clim([_cmin, _cmax])
            img2.set_clim([_cmin, _cmax])
            plt.draw()
        s_cmin.on_changed(update)
        s_cmax.on_changed(update)
        resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
        button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')
        def reset(event):
            s_cmin.reset()
            s_cmax.reset()
        button.on_clicked(reset)
        plt.show()


def open_multiple_fits():
    fileName = filedialog.askopenfilenames(title = "Select file",
    filetypes = (("FIT files","*.fit"),("FITS files","*.fits")))
    calibrated = False
    main_counter = 0

    for file in fileName:
        with open("Values.csv", mode='a+', newline='') as csv_file:
            writer=csv.writer(csv_file)
            writer.writerow([file])
        csv_filenames.append(file)
        opened_file = fits.open(file)
        data = opened_file[0]
        if(main_counter==0):
            calibrated = msb.askquestion("Kablibracja", "Czy wczytane FITSy wymagają kalibracji?")
            if calibrated == 'yes':
                calibration(file,main_counter)
                break
            file_show = fits.open(file)[0]
            img_data = file_show.data
            fig = plt.figure()
            plt.imshow(img_data,cmap='gray',  norm=LogNorm())
            plt.colorbar()
            fig.suptitle(fileName)
            canvas = FigureCanvasTkAgg(fig, master=root)
            canvas.mpl_connect('button_press_event', get_cursor_position)
            toolbar = NavigationToolbar2Tk(canvas, root, pack_toolbar=False)
            toolbar.pack(side=tk.BOTTOM, fill=tk.X)
            toolbar.update()
            canvas.draw()
            canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
            #msb.showinfo("Informacja","Wybierz gwiazdę główną dla której chcesz określić tranzyt")
            var = tk.IntVar()
            var2 = tk.IntVar()
            newWindow2 = Toplevel(root)
            newWindow2.title("Gwiazda główna")
            newWindow2.geometry("200x150+850+200")
            button_1 = tk.Button(newWindow2, text="Wybierz gwiazdę główną", command=lambda: var.set(1))
            button_1.place(x=15,y=20)
            button_1.wait_variable(var)
            r = 0
            R = 0
            r,R = open_new_window()
            button_2 =tk.Button(newWindow2, text="Zatwierdź wybór gwiazdy", command=lambda: var2.set(1))
            button_2.place(x=15,y=60)
            button_2.wait_variable(var2)
            if(var2!=0):
                choice = msb.askquestion("Potwierdzenie wyboru", "Czy potwierdzasz wybór gwiazdy głównej?")
        if choice =='yes':
            main_counter += 1
            temp_temp_1, temp_temp_2 = (ring_estimation(file,r,R,srodek_y,srodek_x,main_counter))
            csv_ring_values.append(ring_estimation(file,r,R,srodek_y,srodek_x,main_counter))
            with open("RING_ESTIMATION.csv", mode='a+') as csv_file:
                writer=csv.writer(csv_file)
                writer.writerow([file])
            with open("RING_ESTIMATION.csv", mode='a+') as csv_file:
                writer=csv.writer(csv_file)
                writer.writerow([ring_estimation(file,r,R,srodek_y,srodek_x,main_counter)])

            main_star_brigtness_1 = float(temp_temp_1)
            tlo = float(temp_temp_2)
            main_star_brightness.append(main_star_brigtness_1)
            csv_main_brightness_values.append(main_star_brigtness_1)
            with open("Main_star_brightness.csv", mode='a+', newline='') as csv_file:
                writer=csv.writer(csv_file)
                writer.writerow([file])
            with open("Main_star_brightness.csv", mode='a+', newline='') as csv_file:
                writer=csv.writer(csv_file)
                writer.writerow([main_star_brigtness_1])
            main_background_brightness.append(tlo)
            csv_background_values.append(tlo)
            with open("Values.csv", mode='a+', newline='') as csv_file:
                writer=csv.writer(csv_file)
                writer.writerow([tlo])
            reference_stars_fixed = []
            reference_stars = reference_star_finder(file,srodek_x,srodek_y)
            data_x = []
            data_y = []
            for i in reference_stars:
                data_x.append(i['xcentroid'])
                data_y.append(i['ycentroid'])
            fixed_data_x = []
            fixed_data_y = []
            counter = 0
            for i in range(0,len(data_x)-1):
                if(abs(data_x[counter]-data_x[counter+1])>6 and abs(data_y[counter]-data_y[counter+1])>6):
                    fixed_data_x.append(int(data_x[counter]))
                    fixed_data_y.append(int(data_y[counter]))
                    csv_ref_coords.append([data_x[counter],data_y[counter]])
                    with open("Values.csv", mode='a+', newline='') as csv_file:
                        writer=csv.writer(csv_file)
                        writer.writerow([data_x[counter],data_y[counter]])
                counter += 1
            temp_counter_1 = 0
            for i in fixed_data_x:
                coord_x = i
                coord_y = fixed_data_y[temp_counter_1]
                if(coord_x>100 and coord_y>100):
                    reference_stars_fixed.append((i,fixed_data_y[temp_counter_1]))
                temp_counter_1 += 1
            reference_stars_brightness = []
            for i in reference_stars_fixed:
                coord_x, coord_y = i
                temp = reference_ring(file,coord_x,coord_y,tlo,main_counter)
                with open("reference_ring_values.csv", mode='a+', newline='') as csv_file:
                    writer=csv.writer(csv_file)
                    writer.writerow([temp])
                csv_ref_ring_values.append(temp)
                with open("Ref_ring_values.csv", mode='a+', newline='') as csv_file:
                    writer=csv.writer(csv_file)
                    writer.writerow([file])
                with open("Ref_ring_values.csv", mode='a+', newline='') as csv_file:
                    writer=csv.writer(csv_file)
                    writer.writerow([temp])
                reference_stars_brightness.append(temp)
            temp_reference_stars_brightness = []
            for i in reference_stars_brightness:
                if i:
                    temp_reference_stars_brightness.append(i)
            niewiemcotoda = sum(temp_reference_stars_brightness)
            main_value = main_star_brigtness_1/niewiemcotoda

            print("FITS: ", main_counter)
            print("SUMA temp_reference_stars_brightness: ", niewiemcotoda)
            print("Main_star_brightness: ", main_star_brigtness_1)
            print("Main_value: ", main_value)
            pomiar_gwiazdy_znormalizowany.append(main_value)
            csv_main_values.append(main_value)
            with open("Values.csv", mode='a+', newline='') as csv_file:
                writer=csv.writer(csv_file)
                writer.writerow([main_value])
            with open("Main_values.csv", mode='a+', newline='') as csv_file:
                writer=csv.writer(csv_file)
                writer.writerow([main_value])
            with open("Main_values.csv", mode='a+', newline='') as csv_file:
                writer=csv.writer(csv_file)
                writer.writerow([file])
            print("------------------------------liczenie: ",main_value)
            opened_file.close()

        else: print("er")
    with open("main_value.csv", mode='w') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow("FILENAMES: ")
        for i in csv_filenames:
            writer.writerow([i])
        writer.writerow("MAIN VALUES: ")
        for i in csv_main_values:
            writer.writerow([round(i,5)])
        writer.writerow("RING VALUES: ")
        for i in csv_ring_values:
            writer.writerow([i])
        writer.writerow("MAIN BRIGHTNESS VALUES: ")
        for i in csv_main_brightness_values:
            writer.writerow([i])
        writer.writerow("BACKGROUND VALUES: ")
        for i in csv_background_values:
            writer.writerow([i])
        writer.writerow("REFERENCE RING VALUES: ")
        for i in csv_ref_ring_values:
            writer.writerow([i])
    newWindow1 = Toplevel(root)
    newWindow1.title("Wykres")
    newWindow1.geometry("799x695+800+200")
    normalized_data = []
    for i in pomiar_gwiazdy_znormalizowany:
        if i:
            temp = (-i)+1
            normalized_data.append(temp)
    fig1 = plt.figure()
    #axx = np.arange(0,len(normalized_data))
    #plt.plot(axx,normalized_data)
    axx = np.arange(0,len(normalized_data))
    plt.plot(axx,normalized_data)
    plt.ylim(0.75,1.1)
    #plt.gca().invert_yaxis()
    fig1.suptitle("Wykres")
    canvas1 = FigureCanvasTkAgg(fig1, master=newWindow1)
    canvas1.get_tk_widget().pack(side = tk.BOTTOM, fill = tk.BOTH, expand = 1)

menubar = Menu(root)
filemenu = Menu(menubar, tearoff=0)
filemenu.add_command(label="Open FITS", command=open_fits)
filemenu.add_separator()
filemenu.add_command(label="Open multiple FITS", command =open_multiple_fits)
filemenu.add_separator()
filemenu.add_command(label="Exit", command=root.quit)
menubar.add_cascade(label="File", menu=filemenu)
root.config(menu=menubar)
root.mainloop()
