/****************************************************************************
 *
 *      imagefilm.cc: image data handling class
 *      This is part of the yafray package
 *		See AUTHORS for more information
 *
 *      This library is free software; you can redistribute it and/or
 *      modify it under the terms of the GNU Lesser General Public
 *      License as published by the Free Software Foundation; either
 *      version 2.1 of the License, or (at your option) any later version.
 *
 *      This library is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *      Lesser General Public License for more details.
 *
 *      You should have received a copy of the GNU Lesser General Public
 *      License along with this library; if not, write to the Free Software
 *      Foundation,Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

#include <core_api/imagefilm.h>
#include <core_api/imagehandler.h>
#include <core_api/scene.h>
#include <yafraycore/monitor.h>
#include <yafraycore/timer.h>
#include <utilities/math_utils.h>
#include <resources/yafLogoTiny.h>

#include <cstring>
#include <string>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <iomanip>
#include <utility>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp> 
#include <boost/filesystem.hpp>

#if HAVE_FREETYPE
#include <resources/guifont.h>
#include <ft2build.h>
#include <utilities/stringUtils.h>
#include FT_FREETYPE_H
#endif

__BEGIN_YAFRAY

#define FILTER_TABLE_SIZE 16
#define MAX_FILTER_SIZE 8

//! Simple alpha blending
#define alphaBlend(b_bg_col, b_fg_col, b_alpha) (( b_bg_col * (1.f - b_alpha) ) + ( b_fg_col * b_alpha ))

typedef float filterFunc(float dx, float dy);

float Box(float dx, float dy){ return 1.f; }

/*!
Mitchell-Netravali constants
with B = 1/3 and C = 1/3 as suggested by the authors
mnX1 = constants for 1 <= |x| < 2
mna1 = (-B - 6 * C)/6
mnb1 = (6 * B + 30 * C)/6
mnc1 = (-12 * B - 48 * C)/6
mnd1 = (8 * B + 24 * C)/6

mnX2 = constants for 1 > |x|
mna2 = (12 - 9 * B - 6 * C)/6
mnb2 = (-18 + 12 * B + 6 * C)/6
mnc2 = (6 - 2 * B)/6

#define mna1 -0.38888889
#define mnb1  2.0
#define mnc1 -3.33333333
#define mnd1  1.77777778

#define mna2  1.16666666
#define mnb2 -2.0
#define mnc2  0.88888889
*/
#define gaussExp 0.00247875

float Mitchell(float dx, float dy)
{
	float x = 2.f * fSqrt(dx*dx + dy*dy);

	if(x >= 2.f) return (0.f);

	if(x >= 1.f) // from mitchell-netravali paper 1 <= |x| < 2
	{
		return (float)( x * ( x * ( x * -0.38888889f + 2.0f) - 3.33333333f) + 1.77777778f );
	}

	return (float)( x * x * ( 1.16666666f * x - 2.0f ) + 0.88888889f );
}

float Gauss(float dx, float dy)
{
	float r2 = dx*dx + dy*dy;
	return std::max(0.f, float(fExp(-6 * r2) - gaussExp));
}

//Lanczos sinc window size 2
float Lanczos2(float dx, float dy)
{
	float x = fSqrt(dx*dx + dy*dy);

	if(x == 0.f) return 1.f;

	if(-2 < x && x < 2)
	{
		float a = M_PI * x;
		float b = M_PI_2 * x;
		return ((fSin(a) * fSin(b)) / (a*b));
	}

	return 0.f;
}

imageFilm_t::imageFilm_t (int width, int height, int xstart, int ystart, colorOutput_t &out, float filterSize, filterType filt,
						  renderEnvironment_t *e, bool showSamMask, int tSize, imageSpliter_t::tilesOrderType tOrder, bool pmA):
	w(width), h(height), cx0(xstart), cy0(ystart), filterw(filterSize*0.5), output(&out),
	env(e), showMask(showSamMask), tileSize(tSize), tilesOrder(tOrder), premultAlpha(pmA)
{
	cx1 = xstart + width;
	cy1 = ystart + height;
	filterTable = new float[FILTER_TABLE_SIZE * FILTER_TABLE_SIZE];

    const renderPasses_t * renderPasses = env->getRenderPasses();
    
	//Creation of the image buffers for the render passes
	for(int idx = 0; idx < renderPasses->extPassesSize(); ++idx)
	{
		imagePasses.push_back(new rgba2DImage_t(width, height));
	}

	//Creation of the image buffers for the auxiliary render passes
	for(int idx = 0; idx < renderPasses->auxPassesSize(); ++idx)
	{
		auxImagePasses.push_back(new rgba2DImage_t(width, height));
	}

	densityImage = nullptr;
	estimateDensity = false;
	dpimage = nullptr;

	// fill filter table:
	float *fTp = filterTable;
	float scale = 1.f/(float)FILTER_TABLE_SIZE;

	filterFunc *ffunc=0;
	switch(filt)
	{
		case MITCHELL: ffunc = Mitchell; filterw *= 2.6f; break;
		case LANCZOS: ffunc = Lanczos2; break;
		case GAUSS: ffunc = Gauss; filterw *= 2.f; break;
		case BOX:
		default:	ffunc = Box;
	}

	filterw = std::min(std::max(0.501f, filterw), 0.5f * MAX_FILTER_SIZE); // filter needs to cover at least the area of one pixel and no more than MAX_FILTER_SIZE/2

	for(int y=0; y < FILTER_TABLE_SIZE; ++y)
	{
		for(int x=0; x < FILTER_TABLE_SIZE; ++x)
		{
			*fTp = ffunc((x+.5f)*scale, (y+.5f)*scale);
			++fTp;
		}
	}

	tableScale = 0.9999 * FILTER_TABLE_SIZE/filterw;
	area_cnt = 0;

	pbar = new ConsoleProgressBar_t(80);
	session.setStatusCurrentPassPercent(pbar->getPercent());
	
	AA_detect_color_noise = false;
	AA_dark_threshold_factor = 0.f;
	AA_variance_edge_size = 10;
	AA_variance_pixels = 0;
	AA_clamp_samples = 0.f;	
}

imageFilm_t::~imageFilm_t ()
{
	
	//Deletion of the image buffers for the additional render passes
	for(size_t idx = 0; idx < imagePasses.size(); ++idx)
	{
		delete(imagePasses[idx]);
	}
	imagePasses.clear();

	for(size_t idx = 0; idx < auxImagePasses.size(); ++idx)
	{
		delete(auxImagePasses[idx]);
	}
	auxImagePasses.clear();
	
	//Deletion of the auxiliary image buffers
	for(size_t idx = 0; idx < auxImagePasses.size(); ++idx)
	{
		delete(auxImagePasses[idx]);
	}
	auxImagePasses.clear();

	if(densityImage) delete densityImage;
	delete[] filterTable;
	if(splitter) delete splitter;
	if(dpimage) delete dpimage;
	if(pbar) delete pbar; //remove when pbar no longer created by imageFilm_t!!
}

void imageFilm_t::init(int numPasses)
{
	// Clear color buffers
	for(size_t idx = 0; idx < imagePasses.size(); ++idx)
	{
		imagePasses[idx]->clear();
	}

	// Clear density image
	if(estimateDensity)
	{
		if(!densityImage) densityImage = new rgb2DImage_nw_t(w, h);
		else densityImage->clear();
	}

	// Setup the bucket splitter
	if(split)
	{
		next_area = 0;
		scene_t *scene = env->getScene();
		int nThreads = 1;
		if(scene) nThreads = scene->getNumThreads();
		splitter = new imageSpliter_t(w, h, cx0, cy0, tileSize, tilesOrder, nThreads);
		area_cnt = splitter->size();
	}
	else area_cnt = 1;

	if(pbar) pbar->init(w * h);
	session.setStatusCurrentPassPercent(pbar->getPercent());

	abort = false;
	completed_cnt = 0;
	nPass = 1;
	nPasses = numPasses;

	imagesAutoSavePassCounter = 0;
	filmAutoSavePassCounter = 0;
	resetImagesAutoSaveTimer();
	resetFilmAutoSaveTimer();
	gTimer.addEvent("imagesAutoSaveTimer");
	gTimer.addEvent("filmAutoSaveTimer");
	gTimer.start("imagesAutoSaveTimer");
	gTimer.start("filmAutoSaveTimer");

	if(!output->isPreview())	// Avoid doing the Film Load & Save operations and updating the film check values when we are just rendering a preview!
	{
		if(filmFileSaveLoad == FILM_FILE_LOAD_SAVE) imageFilmLoadAllInFolder();	//Load all the existing Film in the images output folder, combining them together. It will load only the Film files with the same "base name" as the output image film (including file name, computer node name and frame) to allow adding samples to animations.
	
		if(filmFileSaveLoad == FILM_FILE_LOAD_SAVE || filmFileSaveLoad == FILM_FILE_SAVE) imageFilmFileBackup(); //If the imageFilm is set to Save, at the start rename the previous film file as a "backup" just in case the user has made a mistake and wants to get the previous film back. 

		imageFilmUpdateCheckInfo(); //film load check data initialization. Make sure this is done after the Film Load operation (if any).
	}
}

int imageFilm_t::nextPass(int numView, bool adaptive_AA, std::string integratorName, bool skipNextPass)
{
	splitterMutex.lock();
	next_area = 0;
	splitterMutex.unlock();
	nPass++;
	imagesAutoSavePassCounter++;
	filmAutoSavePassCounter++;
	
	if(skipNextPass) return 0;
	
	std::stringstream passString;

	Y_DEBUG << "nPass=" << nPass << " imagesAutoSavePassCounter="<<imagesAutoSavePassCounter<<" filmAutoSavePassCounter="<<filmAutoSavePassCounter<<yendl;

	if(session.renderInProgress() && !output->isPreview())	//avoid saving images/film if we are just rendering material/world/lights preview windows, etc
	{
		colorOutput_t *out2 = env->getOutput2();

		if((imagesAutoSaveIntervalType == AUTOSAVE_PASS_INTERVAL) && (imagesAutoSavePassCounter >= imagesAutoSaveIntervalPasses))
		{
			if(output && output->isImageOutput()) this->flush(numView, IF_ALL, output);
			else if(out2 && out2->isImageOutput()) this->flush(numView, IF_ALL, out2);
			imagesAutoSavePassCounter = 0;
		}

		if((filmFileSaveLoad == FILM_FILE_LOAD_SAVE || filmFileSaveLoad == FILM_FILE_SAVE) && (filmAutoSaveIntervalType == AUTOSAVE_PASS_INTERVAL) && (filmAutoSavePassCounter >= filmAutoSaveIntervalPasses))
		{
			if((output && output->isImageOutput()) || (out2 && out2->isImageOutput()))
			{
				imageFilmSave();
				filmAutoSavePassCounter = 0;
			}
		}
	}

	const renderPasses_t * renderPasses = env->getRenderPasses();

    rgba2DImage_t * samplingFactorImagePass = getImagePassFromIntPassType(PASS_INT_DEBUG_SAMPLING_FACTOR);
	
    if(flags) flags->clear();
	else flags = new tiledBitArray2D_t<3>(w, h, true);
    std::vector<colorA_t> colExtPasses(imagePasses.size(), colorA_t(0.f));
	int variance_half_edge = AA_variance_edge_size / 2;
	
	float AA_thresh_scaled = AA_thesh;

	int n_resample=0;
	
	if(adaptive_AA && AA_thesh > 0.f)
	{
		for(int y=0; y<h-1; ++y)
		{
			for(int x = 0; x < w-1; ++x)
			{
				flags->clearBit(x, y);
			}
		}
		        
		for(int y=0; y<h-1; ++y)
		{
			for(int x = 0; x < w-1; ++x)
			{
                //We will only consider the Combined Pass (pass 0) for the AA additional sampling calculations.

				if((*imagePasses.at(0))(x, y).weight <= 0.f) flags->setBit(x, y);	//If after reloading ImageFiles there are pixels that were not yet rendered at all, make sure they are marked to be rendered in the next AA pass

                float matSampleFactor = 1.f;
                if(samplingFactorImagePass)
				{
					matSampleFactor = (*samplingFactorImagePass)(x, y).normalized().R;
					if(!backgroundResampling && matSampleFactor == 0.f) continue;
				}

				colorA_t pixCol = (*imagePasses.at(0))(x, y).normalized();
				float pixColBri = pixCol.abscol2bri();

				if(AA_dark_detection_type == DARK_DETECTION_LINEAR && AA_dark_threshold_factor > 0.f)
				{
					if(AA_dark_threshold_factor > 0.f) AA_thresh_scaled = AA_thesh*((1.f-AA_dark_threshold_factor) + (pixColBri*AA_dark_threshold_factor));
				}
				else if(AA_dark_detection_type == DARK_DETECTION_CURVE)
				{
					AA_thresh_scaled = dark_threshold_curve_interpolate(pixColBri);
				}
				
				if(pixCol.colorDifference((*imagePasses.at(0))(x+1, y).normalized(), AA_detect_color_noise) >= AA_thresh_scaled)
				{
					flags->setBit(x, y); flags->setBit(x+1, y);
				}
				if(pixCol.colorDifference((*imagePasses.at(0))(x, y+1).normalized(), AA_detect_color_noise) >= AA_thresh_scaled)
				{
					flags->setBit(x, y); flags->setBit(x, y+1);
				}
				if(pixCol.colorDifference((*imagePasses.at(0))(x+1, y+1).normalized(), AA_detect_color_noise) >= AA_thresh_scaled)
				{
					flags->setBit(x, y); flags->setBit(x+1, y+1);
				}
				if(x > 0 && pixCol.colorDifference((*imagePasses.at(0))(x-1, y+1).normalized(), AA_detect_color_noise) >= AA_thresh_scaled)
				{
					flags->setBit(x, y); flags->setBit(x-1, y+1);
				}
				
				if(AA_variance_pixels > 0)
				{
					int variance_x = 0, variance_y = 0;//, pixelcount = 0;
					
					//float window_accum = 0.f, window_avg = 0.f;
					
					for(int xd = -variance_half_edge; xd < variance_half_edge - 1 ; ++xd)
					{
						int xi = x + xd;
						if(xi<0) xi = 0;
						else if(xi>=w-1) xi = w-2;
						
						colorA_t cx0 = (*imagePasses.at(0))(xi, y).normalized();
						colorA_t cx1 = (*imagePasses.at(0))(xi+1, y).normalized();
						
						if(cx0.colorDifference(cx1, AA_detect_color_noise) >= AA_thresh_scaled) ++variance_x;
					}
					
					for(int yd = -variance_half_edge; yd < variance_half_edge - 1 ; ++yd)
					{
						int yi = y + yd;
						if(yi<0) yi = 0;
						else if(yi>=h-1) yi = h-2;
						
						colorA_t cy0 = (*imagePasses.at(0))(x, yi).normalized();
						colorA_t cy1 = (*imagePasses.at(0))(x, yi+1).normalized();
						
						if(cy0.colorDifference(cy1, AA_detect_color_noise) >= AA_thresh_scaled) ++variance_y;
					}

					if(variance_x + variance_y >= AA_variance_pixels)
					{
						for(int xd = -variance_half_edge; xd < variance_half_edge; ++xd)
						{
							for(int yd = -variance_half_edge; yd < variance_half_edge; ++yd)
							{
								int xi = x + xd;
								if(xi<0) xi = 0;
								else if(xi>=w) xi = w-1;

								int yi = y + yd;
								if(yi<0) yi = 0;
								else if(yi>=h) yi = h-1;

								flags->setBit(xi, yi);
							}
						}
					}
				}
			}
		}

		for(int y=0; y<h; ++y)
		{
			for(int x = 0; x < w; ++x)
			{
				if(flags->getBit(x, y))
				{	
					++n_resample;
												
					if(session.isInteractive() && showMask)
					{
                        float matSampleFactor = 1.f;
                        if(samplingFactorImagePass)
                        {
                            matSampleFactor = (*samplingFactorImagePass)(x, y).normalized().R;
                            if(!backgroundResampling && matSampleFactor == 0.f) continue;
                        }
                        
						for(size_t idx = 0; idx < imagePasses.size(); ++idx)
						{
							color_t pix = (*imagePasses[idx])(x, y).normalized();
							float pixColBri = pix.abscol2bri();

							if(pix.R < pix.G && pix.R < pix.B)
								colExtPasses[idx].set(0.7f, pixColBri, matSampleFactor > 1.f ? 0.7f : pixColBri);
							else
								colExtPasses[idx].set(pixColBri, 0.7f, matSampleFactor > 1.f ? 0.7f : pixColBri);
						}
						output->putPixel(numView, x, y, renderPasses, colExtPasses, false);
					}
				}
			}
		}
	}
	else
	{
		n_resample = h*w;
	}

	if(session.isInteractive())	output->flush(numView, renderPasses);

	if(session.renderResumed()) passString << "Film loaded + ";
	
	passString << "Rendering pass " << nPass << " of " << nPasses << ", resampling " << n_resample << " pixels.";

	Y_INFO << integratorName << ": " << passString.str() << yendl;

	if(pbar)
	{
		pbar->init(w * h);
		session.setStatusCurrentPassPercent(pbar->getPercent());
		pbar->setTag(passString.str().c_str());
	}
	completed_cnt = 0;
	
	return n_resample;
}

bool imageFilm_t::nextArea(int numView, renderArea_t &a)
{
	if(abort) return false;

	int ifilterw = (int) ceil(filterw);

	if(split)
	{
		int n;
		splitterMutex.lock();
		n = next_area++;
		splitterMutex.unlock();

		if(	splitter->getArea(n, a) )
		{
			a.sx0 = a.X + ifilterw;
			a.sx1 = a.X + a.W - ifilterw;
			a.sy0 = a.Y + ifilterw;
			a.sy1 = a.Y + a.H - ifilterw;

			if(session.isInteractive())
			{
				outMutex.lock();
				int end_x = a.X+a.W, end_y = a.Y+a.H;
				output->highliteArea(numView, a.X, a.Y, end_x, end_y);
				outMutex.unlock();
			}

			return true;
		}
	}
	else
	{
		if(area_cnt) return false;
		a.X = cx0;
		a.Y = cy0;
		a.W = w;
		a.H = h;
		a.sx0 = a.X + ifilterw;
		a.sx1 = a.X + a.W - ifilterw;
		a.sy0 = a.Y + ifilterw;
		a.sy1 = a.Y + a.H - ifilterw;
		++area_cnt;
		return true;
	}
	return false;
}

void imageFilm_t::finishArea(int numView, renderArea_t &a)
{
	outMutex.lock();

    const renderPasses_t * renderPasses = env->getRenderPasses();
    
	int end_x = a.X+a.W-cx0, end_y = a.Y+a.H-cy0;

    std::vector<colorA_t> colExtPasses(imagePasses.size(), colorA_t(0.f));

	for(int j=a.Y-cy0; j<end_y; ++j)
	{
		for(int i=a.X-cx0; i<end_x; ++i)
		{
			for(size_t idx = 0; idx < imagePasses.size(); ++idx)
			{
				if(renderPasses->intPassTypeFromExtPassIndex(idx) == PASS_INT_AA_SAMPLES)
				{
					colExtPasses[idx] = (*imagePasses[idx])(i, j).weight;
				}
				else if(renderPasses->intPassTypeFromExtPassIndex(idx) == PASS_INT_OBJ_INDEX_ABS ||
                        renderPasses->intPassTypeFromExtPassIndex(idx) == PASS_INT_OBJ_INDEX_AUTO_ABS ||
                        renderPasses->intPassTypeFromExtPassIndex(idx) == PASS_INT_MAT_INDEX_ABS ||
                        renderPasses->intPassTypeFromExtPassIndex(idx) == PASS_INT_MAT_INDEX_AUTO_ABS
                        )
				{
					colExtPasses[idx] = (*imagePasses[idx])(i, j).normalized();
                    colExtPasses[idx].ceil(); //To correct the antialiasing and ceil the "mixed" values to the upper integer
				}
                else
                {
                    colExtPasses[idx] = (*imagePasses[idx])(i, j).normalized();
                }

				colExtPasses[idx].clampRGB0();
				colExtPasses[idx].ColorSpace_from_linearRGB(colorSpace, gamma);//FIXME DAVID: what passes must be corrected and what do not?
				if(premultAlpha && idx == 0) colExtPasses[idx].alphaPremultiply();

				//To make sure we don't have any weird Alpha values outside the range [0.f, +1.f]
				if(colExtPasses[idx].A < 0.f) colExtPasses[idx].A = 0.f;
				else if(colExtPasses[idx].A > 1.f) colExtPasses[idx].A = 1.f;
			}

			if( !output->putPixel(numView, i, j, renderPasses, colExtPasses) ) abort=true;
		}
	}

	for(size_t idx = 1; idx < imagePasses.size(); ++idx)
	{
		if(renderPasses->intPassTypeFromExtPassIndex(idx) == PASS_INT_DEBUG_FACES_EDGES)
		{
			generateDebugFacesEdges(numView, idx, a.X-cx0, end_x, a.Y-cy0, end_y, true, output);
		}
		
		if(renderPasses->intPassTypeFromExtPassIndex(idx) == PASS_INT_DEBUG_OBJECTS_EDGES || renderPasses->intPassTypeFromExtPassIndex(idx) == PASS_INT_TOON)
		{
			generateToonAndDebugObjectEdges(numView, idx, a.X-cx0, end_x, a.Y-cy0, end_y, true, output);
		}
	}

	if(session.isInteractive()) output->flushArea(numView, a.X, a.Y, end_x+cx0, end_y+cy0, renderPasses);
	
	if(session.renderInProgress() && !output->isPreview())	//avoid saving images/film if we are just rendering material/world/lights preview windows, etc
	{
		gTimer.stop("imagesAutoSaveTimer");
		imagesAutoSaveTimer += gTimer.getTime("imagesAutoSaveTimer");
		if(imagesAutoSaveTimer < 0.f) resetImagesAutoSaveTimer(); //to solve some strange very negative value when using yafaray-xml, race condition somewhere?
		gTimer.start("imagesAutoSaveTimer");

		gTimer.stop("filmAutoSaveTimer");
		filmAutoSaveTimer += gTimer.getTime("filmAutoSaveTimer");
		if(filmAutoSaveTimer < 0.f) resetFilmAutoSaveTimer(); //to solve some strange very negative value when using yafaray-xml, race condition somewhere?
		gTimer.start("filmAutoSaveTimer");

		colorOutput_t *out2 = env->getOutput2();

		if((imagesAutoSaveIntervalType == AUTOSAVE_TIME_INTERVAL) && (imagesAutoSaveTimer > imagesAutoSaveIntervalSeconds))
		{
			Y_DEBUG << "imagesAutoSaveTimer="<<imagesAutoSaveTimer<<yendl;
			if(output && output->isImageOutput()) this->flush(numView, IF_ALL, output);
			else if(out2 && out2->isImageOutput()) this->flush(numView, IF_ALL, out2);
			resetImagesAutoSaveTimer();
		}

		if((filmFileSaveLoad == FILM_FILE_LOAD_SAVE || filmFileSaveLoad == FILM_FILE_SAVE) && (filmAutoSaveIntervalType == AUTOSAVE_TIME_INTERVAL) && (filmAutoSaveTimer > filmAutoSaveIntervalSeconds))
		{
			Y_DEBUG << "filmAutoSaveTimer="<<filmAutoSaveTimer<<yendl;
			if((output && output->isImageOutput()) || (out2 && out2->isImageOutput()))
			{
				imageFilmSave();
			}
			resetFilmAutoSaveTimer();
		}
	}

    if(pbar)
    {
        if(++completed_cnt == area_cnt) pbar->done();
        else pbar->update(a.W * a.H);
        session.setStatusCurrentPassPercent(pbar->getPercent());
    }
    
	outMutex.unlock();
}

void imageFilm_t::flush(int numView, int flags, colorOutput_t *out)
{
    const renderPasses_t * renderPasses = env->getRenderPasses();
    
	if(session.renderFinished())
	{
		outMutex.lock();
		Y_INFO << "imageFilm: Flushing buffer (View number " << numView << ")..." << yendl;
	}

	colorOutput_t *out1 = out ? out : output;
	colorOutput_t *out2 = env->getOutput2();

	if(out1->isPreview()) out2 = nullptr;	//disable secondary file output when rendering a Preview (material preview, etc)
	
	if(out1 == out2) out1 = nullptr;	//if we are already flushing the secondary output (out2) as main output (out1), then disable out1 to avoid duplicated work

	std::string version = session.getYafaRayCoreVersion();

	std::stringstream ssBadge;

	if(!yafLog.getLoggingTitle().empty()) ssBadge << yafLog.getLoggingTitle() << "\n";
	if(!yafLog.getLoggingAuthor().empty() && !yafLog.getLoggingContact().empty()) ssBadge << yafLog.getLoggingAuthor() << " | " << yafLog.getLoggingContact() << "\n";
	else if(!yafLog.getLoggingAuthor().empty() && yafLog.getLoggingContact().empty()) ssBadge << yafLog.getLoggingAuthor() << "\n";
	else if(yafLog.getLoggingAuthor().empty() && !yafLog.getLoggingContact().empty()) ssBadge << yafLog.getLoggingContact() << "\n";
	if(!yafLog.getLoggingComments().empty()) ssBadge << yafLog.getLoggingComments() << "\n";

	ssBadge << "\nYafaRay (" << version << ")" << " " << sysInfoGetOS() << sysInfoGetArchitecture() << sysInfoGetPlatform() << sysInfoGetCompiler(); 

	ssBadge << std::setprecision(2);
	double times = gTimer.getTimeNotStopping("rendert");
	if(session.renderFinished()) times = gTimer.getTime("rendert");
	int timem, timeh;
	gTimer.splitTime(times, &times, &timem, &timeh);
	ssBadge << " | " << w << "x" << h;
	if(session.renderInProgress()) ssBadge << " | " << (session.renderResumed() ? "film loaded + " : "") << "in progress " << std::fixed << std::setprecision(1) << session.currentPassPercent() << "% of pass: " << session.currentPass() << " / " << session.totalPasses();
	else if(session.renderAborted()) ssBadge << " | " << (session.renderResumed() ? "film loaded + " : "") << "stopped at " << std::fixed << std::setprecision(1) << session.currentPassPercent() << "% of pass: " << session.currentPass() << " / " << session.totalPasses();
	else 
	{
		if(session.renderResumed())	ssBadge << " | film loaded + " << session.totalPasses()-1 << " passes";
		else ssBadge << " | " << session.totalPasses() << " passes";
	}
	//if(cx0 != 0) ssBadge << ", xstart=" << cx0; 
	//if(cy0 != 0) ssBadge << ", ystart=" << cy0;
	ssBadge << " | Render time:";
	if (timeh > 0) ssBadge << " " << timeh << "h";
	if (timem > 0) ssBadge << " " << timem << "m";
	ssBadge << " " << times << "s";
	
	times = gTimer.getTimeNotStopping("rendert") + gTimer.getTime("prepass");
	if(session.renderFinished()) times = gTimer.getTime("rendert") + gTimer.getTime("prepass");
	gTimer.splitTime(times, &times, &timem, &timeh);
	ssBadge << " | Total time:";
	if (timeh > 0) ssBadge << " " << timeh << "h";
	if (timem > 0) ssBadge << " " << timem << "m";
	ssBadge << " " << times << "s";
	
	std::stringstream ssLog;
	ssLog << ssBadge.str();
	yafLog.setRenderInfo(ssBadge.str());
	
	if(yafLog.getDrawRenderSettings()) ssBadge << " | " << yafLog.getRenderSettings();
	if(yafLog.getDrawAANoiseSettings()) ssBadge << "\n" << yafLog.getAANoiseSettings();
	if(output && output->isImageOutput()) ssBadge << " " << output->getDenoiseParams();
	else if(out2 && out2->isImageOutput()) ssBadge << " " << out2->getDenoiseParams();

	ssLog << " | " << yafLog.getRenderSettings();
	ssLog << "\n" << yafLog.getAANoiseSettings();
	if(output && output->isImageOutput()) ssLog << " " << output->getDenoiseParams();
	else if(out2 && out2->isImageOutput()) ssLog << " " << out2->getDenoiseParams();

	if(yafLog.getUseParamsBadge())
	{
		if((out1 && out1->isImageOutput()) || (out2 && out2->isImageOutput())) drawRenderSettings(ssBadge);
	}

	if(session.renderFinished())
	{
		Y_PARAMS << "--------------------------------------------------------------------------------" << yendl;
		for (std::string line; std::getline(ssLog, line, '\n');) if(line != "" && line != "\n") Y_PARAMS << line << yendl;
		Y_PARAMS << "--------------------------------------------------------------------------------" << yendl;
	}
	
#ifndef HAVE_FREETYPE
	Y_WARNING << "imageFilm: Compiled without FreeType support." << yendl;
	Y_WARNING << "imageFilm: Text on the parameters badge won't be available." << yendl;
#endif

	float densityFactor = 0.f;

	if(estimateDensity && numDensitySamples > 0) densityFactor = (float) (w * h) / (float) numDensitySamples;

    std::vector<colorA_t> colExtPasses(imagePasses.size(), colorA_t(0.f));

    std::vector<colorA_t> colExtPasses2;	//For secondary file output (when enabled)
	if(out2) colExtPasses2.resize(imagePasses.size(), colorA_t(0.f));

	int outputDisplaceRenderedImageBadgeHeight = 0, out2DisplaceRenderedImageBadgeHeight = 0;
	
	if(out1 && out1->isImageOutput() && yafLog.isParamsBadgeTop()) outputDisplaceRenderedImageBadgeHeight = yafLog.getBadgeHeight();
	
	if(out2 && out2->isImageOutput() && yafLog.isParamsBadgeTop()) out2DisplaceRenderedImageBadgeHeight = yafLog.getBadgeHeight();

	for(int j = 0; j < h; j++)
	{
		for(int i = 0; i < w; i++)
		{
			for(size_t idx = 0; idx < imagePasses.size(); ++idx)
			{
				if(renderPasses->intPassTypeFromExtPassIndex(idx) == PASS_INT_AA_SAMPLES)
				{
					colExtPasses[idx] = (*imagePasses[idx])(i, j).weight;
				}
				else if(renderPasses->intPassTypeFromExtPassIndex(idx) == PASS_INT_OBJ_INDEX_ABS ||
                        renderPasses->intPassTypeFromExtPassIndex(idx) == PASS_INT_OBJ_INDEX_AUTO_ABS ||
                        renderPasses->intPassTypeFromExtPassIndex(idx) == PASS_INT_MAT_INDEX_ABS ||
                        renderPasses->intPassTypeFromExtPassIndex(idx) == PASS_INT_MAT_INDEX_AUTO_ABS
                        )
				{
					colExtPasses[idx] = (*imagePasses[idx])(i, j).normalized();
                    colExtPasses[idx].ceil(); //To correct the antialiasing and ceil the "mixed" values to the upper integer
				}
                else
                {
                    if(flags & IF_IMAGE) colExtPasses[idx] = (*imagePasses[idx])(i, j).normalized();
                    else colExtPasses[idx] = colorA_t(0.f);
                }
								
				if(estimateDensity && (flags & IF_DENSITYIMAGE) && idx == 0 && densityFactor > 0.f) colExtPasses[idx] += colorA_t((*densityImage)(i, j) * densityFactor, 0.f);
                
				colExtPasses[idx].clampRGB0();
				
				if(out2) colExtPasses2[idx] = colExtPasses[idx];
				
				colExtPasses[idx].ColorSpace_from_linearRGB(colorSpace, gamma);//FIXME DAVID: what passes must be corrected and what do not?
				
				if(out2) colExtPasses2[idx].ColorSpace_from_linearRGB(colorSpace2, gamma2);

				if(premultAlpha && idx == 0) 
				{
					colExtPasses[idx].alphaPremultiply();
				}

				if(out2 && premultAlpha2 && idx == 0) 
				{
					colExtPasses2[idx].alphaPremultiply();
				}

				//To make sure we don't have any weird Alpha values outside the range [0.f, +1.f]
				if(colExtPasses[idx].A < 0.f) colExtPasses[idx].A = 0.f;
				else if(colExtPasses[idx].A > 1.f) colExtPasses[idx].A = 1.f;

				if(out2)
				{
					if(colExtPasses2[idx].A < 0.f) colExtPasses2[idx].A = 0.f;
					else if(colExtPasses2[idx].A > 1.f) colExtPasses2[idx].A = 1.f;
				}
			}

			if(out1) out1->putPixel(numView, i, j+outputDisplaceRenderedImageBadgeHeight, renderPasses, colExtPasses);
			if(out2) out2->putPixel(numView, i, j+out2DisplaceRenderedImageBadgeHeight, renderPasses, colExtPasses2);
		}
	}

	for(size_t idx = 1; idx < imagePasses.size(); ++idx)
	{
		if(renderPasses->intPassTypeFromExtPassIndex(idx) == PASS_INT_DEBUG_FACES_EDGES)
		{
			generateDebugFacesEdges(numView, idx, 0, w, 0, h, false, out1, outputDisplaceRenderedImageBadgeHeight, out2, out2DisplaceRenderedImageBadgeHeight);
		}
		
		if(renderPasses->intPassTypeFromExtPassIndex(idx) == PASS_INT_DEBUG_OBJECTS_EDGES || renderPasses->intPassTypeFromExtPassIndex(idx) == PASS_INT_TOON)
		{
			generateToonAndDebugObjectEdges(numView, idx, 0, w, 0, h, false, out1, outputDisplaceRenderedImageBadgeHeight, out2, out2DisplaceRenderedImageBadgeHeight);
		}
	}
	

	if(yafLog.getUseParamsBadge() && dpimage)
	{
		if(out1 && out1->isImageOutput())
		{
			int badgeStartY = h;
			
			if(yafLog.isParamsBadgeTop()) badgeStartY = 0;
			
			for(int j = badgeStartY; j < badgeStartY+dpHeight; j++)
			{
				for(int i = 0; i < w; i++)
				{
					for(size_t idx = 0; idx < imagePasses.size(); ++idx)
					{
						colorA_t &dpcol = (*dpimage)(i, j-badgeStartY);
						colExtPasses[idx] = colorA_t(dpcol, 1.f);
					}
					out1->putPixel(numView, i, j, renderPasses, colExtPasses);
				}
			}
		}

		if(out2 && out2->isImageOutput())
		{
			int badgeStartY = h;
			
			if(yafLog.isParamsBadgeTop()) badgeStartY = 0;
			
			for(int j = badgeStartY; j < badgeStartY+dpHeight; j++)
			{
				for(int i = 0; i < w; i++)
				{
					for(size_t idx = 0; idx < imagePasses.size(); ++idx)
					{
						colorA_t &dpcol = (*dpimage)(i, j-badgeStartY);
						colExtPasses2[idx] = colorA_t(dpcol, 1.f);
					}
					out2->putPixel(numView, i, j, renderPasses, colExtPasses2);
				}
			}
		}
	}

	if(out1 && (session.renderFinished() || out1->isImageOutput())) 
	{
		std::stringstream passString;
		if(out1->isImageOutput()) passString << "Saving image files";
		else passString << "Flushing output";

		Y_INFO << passString.str() << yendl;

		std::string oldTag;

		if(pbar)
		{
			oldTag = pbar->getTag();
			pbar->setTag(passString.str().c_str());
		}

		out1->flush(numView, renderPasses);
		
		if(pbar) pbar->setTag(oldTag);
	}
	
	if(out2 && out2->isImageOutput())
	{
		std::stringstream passString;
		passString << "Saving image files";

		Y_INFO << passString.str() << yendl;

		std::string oldTag;

		if(pbar)
		{
			oldTag = pbar->getTag();
			pbar->setTag(passString.str().c_str());
		}

		out2->flush(numView, renderPasses);

		if(pbar) pbar->setTag(oldTag);
	}

	if(session.renderFinished())
	{
		if(!output->isPreview() && (filmFileSaveLoad == FILM_FILE_LOAD_SAVE || filmFileSaveLoad == FILM_FILE_SAVE))
		{
			if((output && output->isImageOutput()) || (out2 && out2->isImageOutput()))
			{
				imageFilmSave();
			}
		}

		gTimer.stop("imagesAutoSaveTimer");
		gTimer.stop("filmAutoSaveTimer");

		yafLog.clearMemoryLog();
		outMutex.unlock();
		Y_VERBOSE << "imageFilm: Done." << yendl;
	}
}

bool imageFilm_t::doMoreSamples(int x, int y) const
{
	return (AA_thesh>0.f) ? flags->getBit(x-cx0, y-cy0) : true;
}

/* CAUTION! Implemantation of this function needs to be thread safe for samples that
	contribute to pixels outside the area a AND pixels that might get
	contributions from outside area a! (yes, really!) */
void imageFilm_t::addSample(colorPasses_t &colorPasses, int x, int y, float dx, float dy, const renderArea_t *a, int numSample, int AA_pass_number, float inv_AA_max_possible_samples)
{
    const renderPasses_t * renderPasses = env->getRenderPasses();
    
	int dx0, dx1, dy0, dy1, x0, x1, y0, y1;

	// get filter extent and make sure we don't leave image area:

	dx0 = std::max(cx0-x,   Round2Int( (double)dx - filterw));
	dx1 = std::min(cx1-x-1, Round2Int( (double)dx + filterw - 1.0));
	dy0 = std::max(cy0-y,   Round2Int( (double)dy - filterw));
	dy1 = std::min(cy1-y-1, Round2Int( (double)dy + filterw - 1.0));

	// get indizes in filter table
	double x_offs = dx - 0.5;

	int xIndex[MAX_FILTER_SIZE+1], yIndex[MAX_FILTER_SIZE+1];

	for (int i=dx0, n=0; i <= dx1; ++i, ++n)
	{
		double d = std::fabs( (double(i) - x_offs) * tableScale);
		xIndex[n] = Floor2Int(d);
	}

	double y_offs = dy - 0.5;

	for (int i=dy0, n=0; i <= dy1; ++i, ++n)
	{
		double d = std::fabs( (double(i) - y_offs) * tableScale);
		yIndex[n] = Floor2Int(d);
	}

	x0 = x+dx0; x1 = x+dx1;
	y0 = y+dy0; y1 = y+dy1;

	imageMutex.lock();

	for (int j = y0; j <= y1; ++j)
	{
		for (int i = x0; i <= x1; ++i)
		{
			// get filter value at pixel (x,y)
			int offset = yIndex[j-y0]*FILTER_TABLE_SIZE + xIndex[i-x0];
			float filterWt = filterTable[offset];

			// update pixel values with filtered sample contribution
			for(size_t idx = 0; idx < imagePasses.size(); ++idx)
			{
				colorA_t col = colorPasses(renderPasses->intPassTypeFromExtPassIndex(idx));
				
				col.clampProportionalRGB(AA_clamp_samples);

				pixel_t &pixel = (*imagePasses[idx])(i - cx0, j - cy0);

				if(premultAlpha) col.alphaPremultiply();

				if(renderPasses->intPassTypeFromExtPassIndex(idx) == PASS_INT_AA_SAMPLES)
				{
					pixel.weight += inv_AA_max_possible_samples / ((x1-x0+1)*(y1-y0+1));
				}
				else
				{
					pixel.col += (col * filterWt);
					pixel.weight += filterWt;
				}
			}
			for(size_t idx = 0; idx < auxImagePasses.size(); ++idx)
			{
				colorA_t col = colorPasses(renderPasses->intPassTypeFromAuxPassIndex(idx));
				
				col.clampProportionalRGB(AA_clamp_samples);

				pixel_t &pixel = (*auxImagePasses[idx])(i - cx0, j - cy0);

				if(premultAlpha) col.alphaPremultiply();

				if(renderPasses->intPassTypeFromAuxPassIndex(idx) == PASS_INT_AA_SAMPLES)
				{
					pixel.weight += inv_AA_max_possible_samples / ((x1-x0+1)*(y1-y0+1));
				}
				else
				{
					pixel.col += (col * filterWt);
					pixel.weight += filterWt;
				}
			}
		}
	}

	imageMutex.unlock();
}

void imageFilm_t::addDensitySample(const color_t& c, int x, int y, float dx, float dy, const renderArea_t *a)
{
	if(!estimateDensity) return;

	int dx0, dx1, dy0, dy1, x0, x1, y0, y1;

	// get filter extent and make sure we don't leave image area:

	dx0 = std::max(cx0-x,   Round2Int( (double)dx - filterw));
	dx1 = std::min(cx1-x-1, Round2Int( (double)dx + filterw - 1.0));
	dy0 = std::max(cy0-y,   Round2Int( (double)dy - filterw));
	dy1 = std::min(cy1-y-1, Round2Int( (double)dy + filterw - 1.0));


	int xIndex[MAX_FILTER_SIZE+1], yIndex[MAX_FILTER_SIZE+1];

	double x_offs = dx - 0.5;
	for (int i=dx0, n=0; i <= dx1; ++i, ++n)
	{
		double d = std::fabs( (double(i) - x_offs) * tableScale);
		xIndex[n] = Floor2Int(d);
	}

	double y_offs = dy - 0.5;
	for (int i=dy0, n=0; i <= dy1; ++i, ++n)
	{
		float d = fabsf((float)((double(i) - y_offs) * tableScale));
		yIndex[n] = Floor2Int(d);
	}

	x0 = x+dx0; x1 = x+dx1;
	y0 = y+dy0; y1 = y+dy1;

	densityImageMutex.lock();

	for (int j = y0; j <= y1; ++j)
	{
		for (int i = x0; i <= x1; ++i)
		{
			int offset = yIndex[j-y0]*FILTER_TABLE_SIZE + xIndex[i-x0];

			color_t &pixel = (*densityImage)(i - cx0, j - cy0);
			pixel += c * filterTable[offset];
		}
	}

	++numDensitySamples;

	densityImageMutex.unlock();
}

void imageFilm_t::setDensityEstimation(bool enable)
{
	if(enable)
	{
		if(!densityImage) densityImage = new rgb2DImage_nw_t(w, h);
		else densityImage->clear();
	}
	else
	{
		if(densityImage) delete densityImage;
	}
	estimateDensity = enable;
}

void imageFilm_t::setColorSpace(colorSpaces_t color_space, float gammaVal)
{
	colorSpace = color_space;
	gamma = gammaVal;
}

void imageFilm_t::setColorSpace2(colorSpaces_t color_space, float gammaVal)
{
	colorSpace2 = color_space;
	gamma2 = gammaVal;
}

void imageFilm_t::setPremult2(bool premult)
{
	premultAlpha2 = premult;
}

void imageFilm_t::setProgressBar(progressBar_t *pb)
{
	if(pbar) delete pbar;
	pbar = pb;
}

void imageFilm_t::setAANoiseParams(bool detect_color_noise, int dark_detection_type, float dark_threshold_factor, int variance_edge_size, int variance_pixels, float clamp_samples)
{
	AA_detect_color_noise = detect_color_noise;
	AA_dark_detection_type = dark_detection_type;
	AA_dark_threshold_factor = dark_threshold_factor;
	AA_variance_edge_size = variance_edge_size;
	AA_variance_pixels = variance_pixels;
	AA_clamp_samples = clamp_samples;
}

#if HAVE_FREETYPE

void imageFilm_t::drawFontBitmap( FT_Bitmap* bitmap, int x, int y)
{
	int i, j, p, q;
	int x_max = std::min(static_cast<int>(x + bitmap->width), dpimage->getWidth());
	int y_max = std::min(static_cast<int>(y + bitmap->rows), dpimage->getHeight());
	color_t textColor(1.f);

	for ( i = x, p = 0; i < x_max; i++, p++ )
	{
		for ( j = y, q = 0; j < y_max; j++, q++ )
		{
			if ( i >= w || j >= h ) continue;

			int tmpBuf = bitmap->buffer[q * bitmap->width + p];

			if (tmpBuf > 0)
			{
				colorA_t &col = (*dpimage)(std::max(0,i), std::max(0,j));
				float alpha = (float) tmpBuf / 255.0;
				col = colorA_t(alphaBlend((color_t)col, textColor, alpha), col.getA());
			}
		}
	}
}

#endif

void imageFilm_t::drawRenderSettings(std::stringstream & ss)
{
	if(dpimage)
	{
		delete dpimage;
		dpimage = nullptr;
	}

	dpHeight = yafLog.getBadgeHeight();

	dpimage = new rgba2DImage_nw_t(w, dpHeight);
#ifdef HAVE_FREETYPE
	FT_Library library;
	FT_Face face;

	FT_GlyphSlot slot;
	FT_Vector pen; // untransformed origin

	std::string text_utf8 = ss.str();
	std::u32string wtext_utf32 = utf8_to_wutf32(text_utf8);

	// set font size at default dpi
	float fontsize = 12.5f * yafLog.getLoggingFontSizeFactor();

	// initialize library
	if (FT_Init_FreeType( &library ))
	{
		Y_ERROR << "imageFilm: FreeType lib couldn't be initialized!" << yendl;
		return;
	}

	// create face object
	std::string fontPath = yafLog.getLoggingFontPath();
	if(fontPath.empty())
	{
		if (FT_New_Memory_Face( library, (const FT_Byte*)guifont, guifont_size, 0, &face ))
		{
			Y_ERROR << "imageFilm: FreeType couldn't load the default font!" << yendl;
			return;
		}
	}
	else if(FT_New_Face( library, fontPath.c_str(), 0, &face ))
	{
		Y_WARNING << "imageFilm: FreeType couldn't load the font '" << fontPath<< "', loading default font." << yendl;
			
		if (FT_New_Memory_Face( library, (const FT_Byte*)guifont, guifont_size, 0, &face ))
		{
			Y_ERROR << "imageFilm: FreeType couldn't load the default font!" << yendl;
			return;
		}
	}

	FT_Select_Charmap(face , ft_encoding_unicode);
	
	// set character size
	if (FT_Set_Char_Size( face, (FT_F26Dot6)(fontsize * 64.0), 0, 0, 0 ))
	{
		Y_ERROR << "imageFilm: FreeType couldn't set the character size!" << yendl;
		return;
	}

	slot = face->glyph;

	int textOffsetY = -1 * (int) ceil(12 * yafLog.getLoggingFontSizeFactor());
	int textInterlineOffset = (int) ceil(13 * yafLog.getLoggingFontSizeFactor());
#endif
	// offsets
	int textOffsetX = 4;

#ifdef HAVE_FREETYPE
	// The pen position in 26.6 cartesian space coordinates
	pen.x = textOffsetX * 64;
	pen.y = textOffsetY * 64;

	// Draw the text
	for ( size_t n = 0; n < wtext_utf32.size(); n++ )
	{
		// Set Coordinates for the carrige return
		if (wtext_utf32[n] == '\n') {
			pen.x = textOffsetX * 64;
			pen.y -= textInterlineOffset * 64;
			fontsize = 9.5f * yafLog.getLoggingFontSizeFactor();
			if (FT_Set_Char_Size( face, (FT_F26Dot6)(fontsize * 64.0), 0, 0, 0 ))
			{
				Y_ERROR << "imageFilm: FreeType couldn't set the character size!" << yendl;
				return;
			}

			continue;
		}

		// Set transformation
		FT_Set_Transform( face, 0, &pen );

		// Load glyph image into the slot (erase previous one)
		if (FT_Load_Char( face, wtext_utf32[n], FT_LOAD_DEFAULT ))
		{
			Y_ERROR << "imageFilm: FreeType Couldn't load the glyph image for: '" << wtext_utf32[n] << "'!" << yendl;
			continue;
		}

		// Render the glyph into the slot
		FT_Render_Glyph( slot, FT_RENDER_MODE_NORMAL );

		// Now, draw to our target surface (convert position)
		drawFontBitmap( &slot->bitmap, slot->bitmap_left, -slot->bitmap_top);

		// increment pen position
		pen.x += slot->advance.x;
		pen.y += slot->advance.y;
	}

	// Cleanup
	FT_Done_Face    ( face );
	FT_Done_FreeType( library );
#endif

	// Draw logo image
	paraMap_t ihParams;
	
	ihParams["for_output"] = false;
	
	bool logo_loaded = false;
	imageHandler_t *logo = nullptr;

	if(!yafLog.getLoggingCustomIcon().empty())
	{
		std::string iconExtension = yafLog.getLoggingCustomIcon().substr(yafLog.getLoggingCustomIcon().find_last_of(".") + 1);
		std::transform(iconExtension.begin(), iconExtension.end(),iconExtension.begin(), ::tolower);
		
		std::string imageHandlerType = "png";
		if(iconExtension == "jpeg") imageHandlerType = "jpg";
		else imageHandlerType = iconExtension;

		ihParams["type"] = imageHandlerType;
		logo = env->createImageHandler("logoLoader", ihParams, false);
	
		if(logo && logo->loadFromFile(yafLog.getLoggingCustomIcon())) logo_loaded = true;
		else Y_WARNING << "imageFilm: custom params badge icon '" << yafLog.getLoggingCustomIcon() << "' could not be loaded. Using default YafaRay icon." << yendl;
	}

	if(!logo_loaded)
	{
		ihParams["type"] = std::string("png");
		logo = env->createImageHandler("logoLoader", ihParams, false);		
		if(logo && logo->loadFromMemory(yafLogoTiny, yafLogoTiny_size)) logo_loaded = true;
		else Y_WARNING << "imageFilm: default YafaRay params badge icon could not be loaded. No icon will be shown." << yendl;
	}
	
	if(logo_loaded)
	{
		if(logo->getWidth() > 80 || logo->getHeight() > 45) Y_WARNING << "imageFilm: custom params badge logo is quite big (" << logo->getWidth() << " x " << logo->getHeight() << "). It could invade other areas in the badge. Please try to keep logo size smaller than 80 x 45, for example." << yendl;
		int lx, ly;
		int imWidth = std::min(logo->getWidth(), w);
		int imHeight = std::min(logo->getHeight(), dpHeight);

		for ( lx = 0; lx < imWidth; lx++ )
			for ( ly = 0; ly < imHeight; ly++ )
				if(yafLog.isParamsBadgeTop()) (*dpimage)(w-imWidth+lx, ly) = logo->getPixel(lx, ly);
				else (*dpimage)(w-imWidth+lx, dpHeight-imHeight+ly) = logo->getPixel(lx, ly);

		delete logo;
	}

	Y_VERBOSE << "imageFilm: Rendering parameters badge created." << yendl;
}

float imageFilm_t::dark_threshold_curve_interpolate(float pixel_brightness)
{
	if(pixel_brightness <= 0.10f) return 0.0001f;
	else if(pixel_brightness > 0.10f && pixel_brightness <= 0.20f) return (0.0001f + (pixel_brightness - 0.10f) * (0.0010f - 0.0001f) / 0.10f);
	else if(pixel_brightness > 0.20f && pixel_brightness <= 0.30f) return (0.0010f + (pixel_brightness - 0.20f) * (0.0020f - 0.0010f) / 0.10f);
	else if(pixel_brightness > 0.30f && pixel_brightness <= 0.40f) return (0.0020f + (pixel_brightness - 0.30f) * (0.0035f - 0.0020f) / 0.10f);
	else if(pixel_brightness > 0.40f && pixel_brightness <= 0.50f) return (0.0035f + (pixel_brightness - 0.40f) * (0.0055f - 0.0035f) / 0.10f);
	else if(pixel_brightness > 0.50f && pixel_brightness <= 0.60f) return (0.0055f + (pixel_brightness - 0.50f) * (0.0075f - 0.0055f) / 0.10f);
	else if(pixel_brightness > 0.60f && pixel_brightness <= 0.70f) return (0.0075f + (pixel_brightness - 0.60f) * (0.0100f - 0.0075f) / 0.10f);
	else if(pixel_brightness > 0.70f && pixel_brightness <= 0.80f) return (0.0100f + (pixel_brightness - 0.70f) * (0.0150f - 0.0100f) / 0.10f);
	else if(pixel_brightness > 0.80f && pixel_brightness <= 0.90f) return (0.0150f + (pixel_brightness - 0.80f) * (0.0250f - 0.0150f) / 0.10f);
	else if(pixel_brightness > 0.90f && pixel_brightness <= 1.00f) return (0.0250f + (pixel_brightness - 0.90f) * (0.0400f - 0.0250f) / 0.10f);
	else if(pixel_brightness > 1.00f && pixel_brightness <= 1.20f) return (0.0400f + (pixel_brightness - 1.00f) * (0.0800f - 0.0400f) / 0.20f);
	else if(pixel_brightness > 1.20f && pixel_brightness <= 1.40f) return (0.0800f + (pixel_brightness - 1.20f) * (0.0950f - 0.0800f) / 0.20f);
	else if(pixel_brightness > 1.40f && pixel_brightness <= 1.80f) return (0.0950f + (pixel_brightness - 1.40f) * (0.1000f - 0.0950f) / 0.40f);
	else return 0.1000f;
}

std::string imageFilm_t::getFilmPath() const
{
	std::string filmPath = session.getPathImageOutput();
	std::stringstream node;
	node << std::setfill('0') << std::setw(4) << computerNode;
	filmPath += " - node " + node.str();
	filmPath += ".film";
	return filmPath;
}

bool imageFilm_t::imageFilmLoad(const std::string &filename)
{
	bool debugXMLformat = false;	//Enable only for debugging purposes

	try
	{
		std::ifstream ifs(filename, std::fstream::binary);
		char *memblock = new char [1];
		ifs.seekg (0, std::ios::beg);
		ifs.read (memblock, 1);
		bool binaryfile = false;
		if(memblock[0] < '0') binaryfile = true;	//If first character in the film file is not an ASCII number then consider it a binary file
		delete [] memblock;
		ifs.seekg (0, std::ios::beg);
		
		if(debugXMLformat)
		{
			Y_INFO << "imageFilm: Loading film from: \"" << filename << "\" in XML format" << yendl;
			boost::archive::xml_iarchive ia(ifs);
			ia >> BOOST_SERIALIZATION_NVP(*this);
			ifs.close();
		}
		else if(binaryfile == true)
		{
			Y_INFO << "imageFilm: Loading film from: \"" << filename << "\" in Binary (non portable) format" << yendl;
			boost::archive::binary_iarchive ia(ifs);
			ia >> BOOST_SERIALIZATION_NVP(*this);
			ifs.close();
		}
		else
		{
			Y_INFO << "imageFilm: Loading film from: \"" << filename << "\" in Text format" << yendl;
			boost::archive::text_iarchive ia(ifs);
			ia >> BOOST_SERIALIZATION_NVP(*this);
			ifs.close();
		}
		Y_VERBOSE << "imageFilm: Film loaded from file." << yendl;
		return true;
	}
	catch(std::exception& ex){
        Y_WARNING << "imageFilm: error '" << ex.what() << "' while loading ImageFilm file: '" << filename << "'" << yendl;
		return false;
    }
}

void imageFilm_t::imageFilmLoadAllInFolder()
{
	std::stringstream passString;
	passString << "Loading ImageFilm files";

	Y_INFO << passString.str() << yendl;

	std::string oldTag;

	if(pbar)
	{
		oldTag = pbar->getTag();
		pbar->setTag(passString.str().c_str());
	}
	
	std::string filmPath = getFilmPath();
	
	std::string baseImageFileName = boost::filesystem::path(session.getPathImageOutput()).stem().string();

	std::string parentPath = boost::filesystem::path(session.getPathImageOutput()).parent_path().string();
	if(parentPath.empty()) parentPath = ".";	//If parent path is empty, set the path to the current folder
	const std::string target_path( parentPath );
	std::vector<std::string> filmFilesList;
	
	try
	{
		boost::filesystem::directory_iterator it_end;
		for(boost::filesystem::directory_iterator it( target_path ); it != it_end; ++it)
		{
			if(!boost::filesystem::is_regular_file(it->status())) continue;
			if(it->path().extension().string() != ".film") continue;
			if(it->path().stem().string().size() < baseImageFileName.size()) continue;
			if(it->path().stem().string().compare(0, baseImageFileName.size(), baseImageFileName) != 0) continue;
			filmFilesList.push_back(it->path().string());
		}
		std::sort(filmFilesList.begin(), filmFilesList.end());

		for(auto filmFile: filmFilesList)
		{
			imageFilm_t *loadedFilm = new imageFilm_t(w, h, cx0, cy0, *output, 1.0, BOX, env);
			loadedFilm->imageFilmLoad(filmFile);
			
			for(size_t idx=0; idx<imagePasses.size(); ++idx)
			{
				for(int i=0; i<w; ++i)
				{
					for(int j=0; j<h; ++j)
					{
						rgba2DImage_t *loadedImageBuffer = loadedFilm->imagePasses[idx];
						(*imagePasses[idx])(i,j).col += (*loadedImageBuffer)(i,j).col;
						(*imagePasses[idx])(i,j).weight += (*loadedImageBuffer)(i,j).weight;
					}
				}
			}
			
			for(size_t idx=0; idx<auxImagePasses.size(); ++idx)
			{
				for(int i=0; i<w; ++i)
				{
					for(int j=0; j<h; ++j)
					{
						rgba2DImage_t *loadedImageBuffer = loadedFilm->auxImagePasses[idx];
						(*auxImagePasses[idx])(i,j).col += (*loadedImageBuffer)(i,j).col;
						(*auxImagePasses[idx])(i,j).weight += (*loadedImageBuffer)(i,j).weight;
					}
				}
			}
			
			if(samplingOffset < loadedFilm->samplingOffset) samplingOffset = loadedFilm->samplingOffset;
			if(baseSamplingOffset < loadedFilm->baseSamplingOffset) baseSamplingOffset = loadedFilm->baseSamplingOffset;
			
			delete loadedFilm;
		}
	}
	catch(const boost::filesystem::filesystem_error& e)
	{
		Y_WARNING << "imageFilm: error during imageFilm loading process: \"" << e.what() << "\"" << yendl;
	}
	
	if(pbar) pbar->setTag(oldTag);
}


bool imageFilm_t::imageFilmSave()
{
	bool debugXMLformat = false;	//Enable only for debugging purposes
	
	std::stringstream passString;
	passString << "Saving internal ImageFilm file";

	Y_INFO << passString.str() << yendl;

	std::string oldTag;

	if(pbar)
	{
		oldTag = pbar->getTag();
		pbar->setTag(passString.str().c_str());
	}

	std::string filmPath = getFilmPath();

	try
	{
		std::ofstream ofs(filmPath+".tmp", std::fstream::binary);

		if(debugXMLformat)
		{
			Y_INFO << "imageFilm: Saving film to: \"" << filmPath << "\" in XML format" << yendl;
			boost::archive::xml_oarchive oa(ofs);
			oa << BOOST_SERIALIZATION_NVP(*this);
			ofs.close();
		}
		else if(filmFileSaveBinaryFormat)
		{
			Y_INFO << "imageFilm: Saving film to: \"" << filmPath << "\" in Binary (non portable) format" << yendl;
			boost::archive::binary_oarchive oa(ofs);
			oa << BOOST_SERIALIZATION_NVP(*this);
			ofs.close();
		}
		else
		{
			Y_INFO << "imageFilm: Saving film to: \"" << filmPath << "\" in Text format" << yendl;
			boost::archive::text_oarchive oa(ofs);
			oa << BOOST_SERIALIZATION_NVP(*this);
			ofs.close();
		}
	Y_VERBOSE << "imageFilm: Film saved to file." << yendl;
	}
	catch(std::exception& ex){
        Y_WARNING << "imageFilm: error '" << ex.what() << "' while saving ImageFilm file: '" << filmPath << "'" << yendl;
		if(pbar) pbar->setTag(oldTag);
		return false;
    }
    
	try
	{
		boost::filesystem::copy_file(filmPath+".tmp", filmPath, boost::filesystem::copy_option::overwrite_if_exists);
		boost::filesystem::remove(filmPath+".tmp");
	}
	catch(const boost::filesystem::filesystem_error& e)
	{
		Y_WARNING << "imageFilm: file operation error \"" << e.what() << yendl;
	}

	if(pbar) pbar->setTag(oldTag);
	
	return true;
}

void imageFilm_t::imageFilmFileBackup() const
{
	std::stringstream passString;
	passString << "Creating backup of the previous ImageFilm file...";

	Y_INFO << passString.str() << yendl;

	std::string oldTag;

	if(pbar)
	{
		oldTag = pbar->getTag();
		pbar->setTag(passString.str().c_str());
	}

	std::string filmPath = getFilmPath();
	std::string filmPathBackup = filmPath+"-previous.bak";
	
	if(boost::filesystem::exists(filmPath))
	{
		Y_VERBOSE << "imageFilm: Creating backup of previously saved film to: \"" << filmPathBackup << "\"" << yendl;
		try
		{	
			boost::filesystem::rename(filmPath, filmPathBackup);
		}
		catch(const boost::filesystem::filesystem_error& e)
		{
			Y_WARNING << "imageFilm: error during imageFilm file backup \"" << e.what() << "\"" << yendl;
		}
	}
	
	if(pbar) pbar->setTag(oldTag);
}

void imageFilm_t::imageFilmUpdateCheckInfo()
{
	filmload_check.filmStructureVersion = FILM_STRUCTURE_VERSION;
	filmload_check.w = w;
	filmload_check.h = h;
	filmload_check.cx0 = cx0;
	filmload_check.cx1 = cx1;
	filmload_check.cy0 = cy0;
	filmload_check.cy1 = cy1;
	filmload_check.numPasses = imagePasses.size();
}

bool imageFilm_t::imageFilmLoadCheckOk() const
{
    const renderPasses_t * renderPasses = env->getRenderPasses();
    
	bool checksOK = true;
	
	if(filmload_check.filmStructureVersion != FILM_STRUCTURE_VERSION)
	{
		checksOK = false;
		Y_WARNING << "imageFilm: loading/reusing film check failed. Film structure version, expected=" << FILM_STRUCTURE_VERSION << ", in reused/loaded film=" << filmload_check.filmStructureVersion << yendl;
	}	
	if(filmload_check.w != w)
	{
		checksOK = false;
		Y_WARNING << "imageFilm: loading/reusing film check failed. Image width, expected=" << w << ", in reused/loaded film=" << filmload_check.w << yendl;
	}
	if(filmload_check.h != h)
	{
		checksOK = false;
		Y_WARNING << "imageFilm: loading/reusing film check failed. Image height, expected=" << h << ", in reused/loaded film=" << filmload_check.h << yendl;
	}
	if(filmload_check.cx0 != cx0)
	{
		checksOK = false;
		Y_WARNING << "imageFilm: loading/reusing film check failed. Border cx0, expected=" << cx0 << ", in reused/loaded film=" << filmload_check.cx0 << yendl;
	}
	if(filmload_check.cx1 != cx1)
	{
		checksOK = false;
		Y_WARNING << "imageFilm: loading/reusing film check failed. Border cx1, expected=" << cx1 << ", in reused/loaded film=" << filmload_check.cx1 << yendl;
	}
	if(filmload_check.cy0 != cy0)
	{
		checksOK = false;
		Y_WARNING << "imageFilm: loading/reusing film check failed. Border cy0, expected=" << cy0 << ", in reused/loaded film=" << filmload_check.cy0 << yendl;
	}
	if(filmload_check.cy1 != cy1)
	{
		checksOK = false;
		Y_WARNING << "imageFilm: loading/reusing film check failed. Border cy1, expected=" << cy1 << ", in reused/loaded film=" << filmload_check.cy1 << yendl;
	}
	if(filmload_check.numPasses != (size_t) renderPasses->extPassesSize())
	{
		checksOK = false;
		Y_WARNING << "imageFilm: loading/reusing film check failed. Number of render passes, expected=" << renderPasses->extPassesSize() << ", in reused/loaded film=" << filmload_check.numPasses << yendl;
	}
	
	if(!checksOK) Y_WARNING << "imageFilm: loading/reusing film failed because parameters are different. The film will be re-generated." << yendl;
	
	Y_DEBUG << "imageFilm: loading/reusing film check results=" << checksOK << ". Expected: film structure version=" << FILM_STRUCTURE_VERSION << ",w="<<w<<",h="<<h<<",cx="<<cx0<<",cy0="<<cy0<<",cx1="<<cx1<<",cy1="<<cy1<<",numPasses="<<renderPasses->extPassesSize()<<" .In Image File: film structure version="<<filmload_check.filmStructureVersion<<",w="<<filmload_check.w<<",h="<<filmload_check.h<<",cx="<<filmload_check.cx0<<",cy0="<<filmload_check.cy0<<",cx1="<<filmload_check.cx1<<",cy1="<<filmload_check.cy1<<",numPasses="<<filmload_check.numPasses << yendl;
	
	return checksOK;
}


//The next edge detection, debug faces/object edges and toon functions will only work if YafaRay is built with OpenCV support
#ifdef HAVE_OPENCV

void imageFilm_t::edgeImageDetection(std::vector<cv::Mat> & imageMat, float edge_threshold, int edge_thickness, float smoothness) const
{	
	//The result of the edges detection will be stored in the first component image of the vector 
	
	//Calculate edges for the different component images and combine them in the first component image
	for(auto it = imageMat.begin(); it != imageMat.end(); it++)
	{
		cv::Laplacian(*it, *it, -1, 3, 1, 0, cv::BORDER_DEFAULT);
		if(it != imageMat.begin()) imageMat.at(0) = cv::max(imageMat.at(0), *it);
	}
	
	//Get a pure black/white image of the edges
	cv::threshold( imageMat.at(0), imageMat.at(0), edge_threshold, 1.0, 0);

	//Make the edges wider if needed
	if(edge_thickness > 1)
	{
		cv::Mat kernel = cv::Mat::ones( edge_thickness, edge_thickness, CV_32F )/ (float)(edge_thickness*edge_thickness);
		cv::filter2D( imageMat.at(0), imageMat.at(0), -1, kernel );
		cv::threshold( imageMat.at(0), imageMat.at(0), 0.1, 1.0, 0);
	}

	//Soften the edges if needed
	if(smoothness > 0.f) cv::GaussianBlur( imageMat.at(0), imageMat.at(0), cv::Size(3,3), /*sigmaX=*/ smoothness );
}

void imageFilm_t::generateDebugFacesEdges(int numView, int idxPass, int xstart, int width, int ystart, int height, bool drawborder, colorOutput_t * out1, int out1displacement, colorOutput_t * out2, int out2displacement)
{
	const renderPasses_t * renderPasses = env->getRenderPasses();
	const int facesEdgeThickness = renderPasses->facesEdgeThickness;
	const float facesEdgeThreshold = renderPasses->facesEdgeThreshold;
	const float facesEdgeSmoothness = renderPasses->facesEdgeSmoothness;

    rgba2DImage_t * normalImagePass = getImagePassFromIntPassType(PASS_INT_NORMAL_GEOM);
    rgba2DImage_t * zDepthImagePass = getImagePassFromIntPassType(PASS_INT_Z_DEPTH_NORM);

	if(normalImagePass && zDepthImagePass)
	{
		std::vector<cv::Mat> imageMat;
		for(int i=0; i<4; ++i) imageMat.push_back(cv::Mat(h, w, CV_32FC1));

		for(int j=ystart; j<height; ++j)
		{
			for(int i=xstart; i<width; ++i)
			{
				color_t colNormal = (*normalImagePass)(i, j).normalized();
				float zDepth = (*zDepthImagePass)(i, j).normalized().A;

				imageMat.at(0).at<float>(j, i) = colNormal.getR();
				imageMat.at(1).at<float>(j, i) = colNormal.getG();
				imageMat.at(2).at<float>(j, i) = colNormal.getB();
				imageMat.at(3).at<float>(j, i) = zDepth;
			}
		}
		
		edgeImageDetection(imageMat, facesEdgeThreshold, facesEdgeThickness, facesEdgeSmoothness);

		for(int j=ystart; j<height; ++j)
		{
			for(int i=xstart; i<width; ++i)
			{
				color_t colEdge = color_t(imageMat.at(0).at<float>(j, i));

				if(drawborder && (i <= xstart+1 || j <= ystart+1 || i >= width-1-1 || j >= height-1-1)) colEdge = colorA_t(0.5f,0.f,0.f,1.f);

				if(out1) out1->putPixel(numView, i, j+out1displacement, renderPasses, idxPass, colEdge);
				if(out2) out2->putPixel(numView, i, j+out2displacement, renderPasses, idxPass, colEdge);
			}
		}
	}
}

void imageFilm_t::generateToonAndDebugObjectEdges(int numView, int idxPass, int xstart, int width, int ystart, int height, bool drawborder, colorOutput_t * out1, int out1displacement, colorOutput_t * out2, int out2displacement)
{
	const renderPasses_t * renderPasses = env->getRenderPasses();
	const float toonPreSmooth = renderPasses->toonPreSmooth;
	const float toonQuantization = renderPasses->toonQuantization;
	const float toonPostSmooth = renderPasses->toonPostSmooth;
	const color_t toonEdgeColor = color_t(renderPasses->toonEdgeColor[0], renderPasses->toonEdgeColor[1], renderPasses->toonEdgeColor[2]);
	const int objectEdgeThickness = renderPasses->objectEdgeThickness;
	const float objectEdgeThreshold = renderPasses->objectEdgeThreshold;
	const float objectEdgeSmoothness = renderPasses->objectEdgeSmoothness;
	
    rgba2DImage_t * normalImagePass = getImagePassFromIntPassType(PASS_INT_NORMAL_SMOOTH);
    rgba2DImage_t * zDepthImagePass = getImagePassFromIntPassType(PASS_INT_Z_DEPTH_NORM);
	
	if(normalImagePass && zDepthImagePass)
	{
		cv::Mat_<cv::Vec3f> imageMatCombinedVec(h, w, CV_32FC3);
		std::vector<cv::Mat> imageMat;
		for(int i=0; i<4; ++i) imageMat.push_back(cv::Mat(h, w, CV_32FC1));
		color_t colNormal(0.f);
		float zDepth = 0.f;

		for(int j=ystart; j<height; ++j)
		{
			for(int i=xstart; i<width; ++i)
			{
				colNormal = (*normalImagePass)(i, j).normalized();
				zDepth = (*zDepthImagePass)(i, j).normalized().A;

				imageMatCombinedVec(j, i)[0] = (*imagePasses[0])(i, j).normalized().B;
				imageMatCombinedVec(j, i)[1] = (*imagePasses[0])(i, j).normalized().G;
				imageMatCombinedVec(j, i)[2] = (*imagePasses[0])(i, j).normalized().R;
										
				imageMat.at(0).at<float>(j, i) = colNormal.getR();
				imageMat.at(1).at<float>(j, i) = colNormal.getG();
				imageMat.at(2).at<float>(j, i) = colNormal.getB();
				imageMat.at(3).at<float>(j, i) = zDepth;
			}
		}

		cv::GaussianBlur( imageMatCombinedVec, imageMatCombinedVec, cv::Size(3,3), toonPreSmooth );
		
		if(toonQuantization > 0.f)
		{
			for(int j=ystart; j<height; ++j)
			{
				for(int i=xstart; i<width; ++i)
				{
					{
						float h, s, v;
						color_t col(imageMatCombinedVec(j, i)[0], imageMatCombinedVec(j, i)[1], imageMatCombinedVec(j, i)[2]);
						col.rgb_to_hsv(h, s, v);
						h = std::round(h / toonQuantization) * toonQuantization;
						s = std::round(s / toonQuantization) * toonQuantization;
						v = std::round(v / toonQuantization) * toonQuantization;
						col.hsv_to_rgb(h, s, v);
						imageMatCombinedVec(j, i)[0] = col.R;
						imageMatCombinedVec(j, i)[1] = col.G;
						imageMatCombinedVec(j, i)[2] = col.B;
					}
				}
			}
			cv::GaussianBlur( imageMatCombinedVec, imageMatCombinedVec, cv::Size(3,3), toonPostSmooth );
		}

		edgeImageDetection(imageMat, objectEdgeThreshold, objectEdgeThickness, objectEdgeSmoothness);

        int idxToon = getImagePassIndexFromIntPassType(PASS_INT_TOON);

		for(int j=ystart; j<height; ++j)
		{
			for(int i=xstart; i<width; ++i)
			{
				float edgeValue = imageMat.at(0).at<float>(j, i);
				color_t colEdge = color_t(edgeValue);

				if(drawborder && (i <= xstart+1 || j <= ystart+1 || i >= width-1-1 || j >= height-1-1)) colEdge = colorA_t(0.5f,0.f,0.f,1.f);
				
				if(out1) out1->putPixel(numView, i, j+out1displacement, renderPasses, idxPass, colEdge);
				if(out2) out2->putPixel(numView, i, j+out2displacement, renderPasses, idxPass, colEdge);
				
				color_t colToon(imageMatCombinedVec(j, i)[2], imageMatCombinedVec(j, i)[1], imageMatCombinedVec(j, i)[0]);
				colToon.blend(toonEdgeColor, edgeValue);
				
				if(drawborder && (i <= xstart+1 || j <= ystart+1 || i >= width-1-1 || j >= height-1-1)) colToon = colorA_t(0.5f,0.f,0.f,1.f);
				
				if(idxToon != -1)
				{
					if(out1)
					{
						colToon.ColorSpace_from_linearRGB(colorSpace, gamma);
						out1->putPixel(numView, i, j+out1displacement, renderPasses, idxToon, colToon);
					}
					if(out2)
					{
						colToon.ColorSpace_from_linearRGB(colorSpace2, gamma2);
						out2->putPixel(numView, i, j+out2displacement, renderPasses, idxToon, colToon);
					}
				}
			}
		}
	}
}

#else   //If not built with OpenCV, these functions will do nothing

void imageFilm_t::generateToonAndDebugObjectEdges(int numView, int idxPass, int xstart, int width, int ystart, int height, bool drawborder, colorOutput_t * out1, int out1displacement, colorOutput_t * out2, int out2displacement) { }

void imageFilm_t::generateDebugFacesEdges(int numView, int idxPass, int xstart, int width, int ystart, int height, bool drawborder, colorOutput_t * out1, int out1displacement, colorOutput_t * out2, int out2displacement) { }

#endif


rgba2DImage_t * imageFilm_t::getImagePassFromIntPassType(int intPassType)
{
    for(size_t idx = 1; idx < imagePasses.size(); ++idx)
	{
		if(env->getScene()->getRenderPasses()->intPassTypeFromExtPassIndex(idx) == intPassType) return imagePasses[idx];
	}
    
	for(size_t idx = 0; idx < auxImagePasses.size(); ++idx)
	{
		if(env->getScene()->getRenderPasses()->intPassTypeFromAuxPassIndex(idx) == intPassType) return auxImagePasses[idx];
	}
    
    return nullptr;
}

int imageFilm_t::getImagePassIndexFromIntPassType(int intPassType)
{
    for(size_t idx = 1; idx < imagePasses.size(); ++idx)
	{
		if(env->getScene()->getRenderPasses()->intPassTypeFromExtPassIndex(idx) == intPassType) return (int) idx;
	}

    return -1;
}

int imageFilm_t::getAuxImagePassIndexFromIntPassType(int intPassType)
{
	for(size_t idx = 0; idx < auxImagePasses.size(); ++idx)
	{
		if(env->getScene()->getRenderPasses()->intPassTypeFromAuxPassIndex(idx) == intPassType) return (int) idx;
	}

    return -1;
}


__END_YAFRAY
