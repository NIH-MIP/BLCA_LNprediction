'''
this inferencer model is now flexible to accept classification
predictions from binary or multi-class predictions
it produces a csv with per-image predictions and softmax outputs
for every class, as well as overall class prediction
'''


from fastai.vision import *
import datetime

class Inferencer:

    def __init__(self):
        self.imagedir = '/media/mip/14cd9269-81aa-49d1-b074-939a9565b962/Stephanie_Harmon/bladder_path/data/decon/20x'
        self.outdir = '/media/mip/14cd9269-81aa-49d1-b074-939a9565b962/Stephanie_Harmon/bladder_path/data/decon/20x/training_log'
        self.testPath = os.path.join(self.imagedir, 'test_NIH')
        self.model_name = '20x_obsv_one5loss_09132019-1431.pkl'
        self.device=1

    def convert_model_to_export(self):
        initial_filename = os.path.join(self.outdir, 'exported_models', self.model_name)
        final_filename = os.path.join(self.outdir, 'exported_models', 'export.pkl')
        shutil.copy2(initial_filename, final_filename)


    def fastai_apply_model(self):
        '''
        applies model on patient level
        :param index(int): where to index in
        :return:
        '''

        torch.cuda.set_device(self.device)

        # set up the output directory
        if not os.path.isdir(self.outdir):
            os.mkdir(self.outdir)
        if not os.path.isdir(os.path.join(self.outdir, 'val_results')):
            os.mkdir(os.path.join(self.outdir, 'val_results'))

        # change this part
        model_path = os.path.join(self.outdir, 'exported_models')
        learn = load_learner(model_path)
        test_path = self.testPath

        df_out = pd.DataFrame()

        for patient in os.listdir(os.path.join(test_path)):
            print(patient)
            for image in os.listdir(os.path.join(test_path, patient)):
                img = open_image(os.path.join(test_path, patient, image))
                pred_class, pred_idx, outputs = learn.predict(img)
                outputs_np = outputs.numpy()
                all_classes = learn.data.classes
                #t_df = pd.DataFrame([image, pred_class, outputs_np]).transpose()
                df_i = pd.DataFrame([image, pred_class]).transpose()
                for classi in range(0,len(all_classes)):
                    p_df = pd.DataFrame([outputs_np[classi]]).transpose()
                    df_i=pd.concat([df_i,p_df], axis=1)
                df_out=pd.concat([df_out,df_i],axis=0)

        #write dataframe out to csv
        df_out.columns=['img_name','class_pred']+all_classes
        df_out.to_csv(os.path.join(self.outdir, 'val_results', 'val_results_by_img' + '_' + self.model_name + '_' + str(
                datetime.datetime.now().strftime("%m%d%Y-%H%M")) + '.csv'))

if __name__ == '__main__':
    d = Inferencer()
    d.convert_model_to_export()
    d.fastai_apply_model()
