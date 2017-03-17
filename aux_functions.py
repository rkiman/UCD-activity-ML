import matplotlib.pyplot as plt
import math as math	
def plot_full(x,y,title='',xlabel = 'x',ylabel = 'y',textsize = '15',graphtype = ''):
	
	plt.plot(x,y,graphtype)		
	if(title != ''):
		plt.title(title,fontsize=textsize)

	plt.xlabel(xlabel,fontsize=textsize)
	plt.ylabel(ylabel,fontsize=textsize)	
	axes = plt.gca()
	#axes.set_xlim([0.1,1])
	axes.set_ylim([-1,10])
	plt.show()
	
def plot_err(x,y,xer='',yer='',title='',xlabel = 'x',ylabel = 'y',textsize = '15',graphtype = ''):
	
	if(xer == ''):
		plt.errorbar(x,y,yerr=yer,linestyle='None',marker = graphtype)	
	elif(yer == ''):
		plt.errorbar(x,y,xerr=xer,linestyle='None',marker = graphtype)				
	else:
		plt.errorbar(x,y,xerr=xer,yerr=yer,linestyle='None',marker = graphtype)
	if(title != ''):
		plt.title(title,fontsize=textsize)

	plt.xlabel(xlabel,fontsize=textsize)
	plt.ylabel(ylabel,fontsize=textsize)	
	axes = plt.gca()
	axes.set_xlim([0.1,1])
	axes.set_ylim([-1,10])
	plt.show()

def plot_hist(x,nbins,title='',textsize = '15'):

	num_stars = len(x)
	text = r'Number of Stars:' + str(num_stars)
	NBINS = nbins
	energy_hist = plt.hist(x, NBINS)
	plt.yscale('log')
	axes = plt.gca()
	axes.set_ylim([10,13000])
	plt.title(title)
	plt.figtext(0.2,0.85,text,fontsize = textsize)
	plt.show()

def error(Fw,Fs,sigw,sigs):
	error = math.sqrt((1/Fs)**2*sigw**2+(Fw/Fs**2)**2*sigs**2)
	return error
