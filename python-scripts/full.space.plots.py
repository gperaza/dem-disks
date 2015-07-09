import plots_3disks_lib as p3d
import matplotlib.pyplot as plt

def main():
    #p3d.plot_angle_time(plotRelax=True)
    #p3d.plot_coord_time("w", iTime=2, fTime=2.5, show=True)
    #p3d.plot_coord_time("y", iTime=9.9, fTime=10)
    #p3d.hist_links()
    #p3d.plot_supperpose(coord1="y", diskId1="0",
    #                    coord2="y", diskId2="1",
    #                    show=False, iTime=9.9, fTime=10)
    #p3d.link_stats_coord(coord="w", show=True, iTime=9.9, fTime=9.9250)

    iTime = 9
    fTime = 10
    p3d.plot_compare_y(show=False, iTime=iTime, fTime=fTime)
    p3d.plot_center_path(iTime=iTime, fTime=fTime, show=False)
    #p3d.plot_compare_supperpose(show=False, iTime=iTime, fTime=fTime)
    #p3d.plot_compare_links(show=False, iTime=iTime, fTime=fTime)

if __name__ == "__main__":
    main()
